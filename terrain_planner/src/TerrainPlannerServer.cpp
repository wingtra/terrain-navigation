#include <crow.h>  // Crow single-header or multi-file, see https://github.com/CrowCpp/Crow
#include <ompl/base/StateSpace.h>
#include <ompl/geometric/PathGeometric.h>
#include <memory>
#include <string>
#include <fstream>
#include <sstream>

#include "terrain_planner/terrain_ompl_rrt.h"
#include "terrain_planner/planner.h"
#include "terrain_navigation/terrain_map.h"

// A helper to convert a planned path (list of 3D states + yaw) to a JSON array
static crow::json::wvalue pathToJSON(const std::vector<Eigen::Vector3d>& path_points,
                                     const std::vector<double>& path_yaws)
{
    // Accumulate each point into a vector of crow::json::wvalue
    std::vector<crow::json::wvalue> arr;
    arr.reserve(path_points.size());

    for (size_t i = 0; i < path_points.size(); i++)
    {
        crow::json::wvalue pt;
        pt["x"]   = path_points[i].x();
        pt["y"]   = path_points[i].y();
        pt["z"]   = path_points[i].z();
        pt["yaw"] = (i < path_yaws.size()) ? path_yaws[i] : 0.0;
        arr.push_back(std::move(pt));
    }
    // Construct and return a crow::json::wvalue from that vector
    return crow::json::wvalue(std::move(arr));
}

// Utility for computing the yaw from consecutive points
inline double computeYaw(const Eigen::Vector3d& from, const Eigen::Vector3d& to)
{
    double dx = to.x() - from.x();
    double dy = to.y() - from.y();
    return std::atan2(dy, dx);
}

int main(int argc, char* argv[])
{
    // Create a crow::SimpleApp
    crow::SimpleApp app;

    // Define a POST endpoint: /plan
    CROW_ROUTE(app, "/plan").methods(crow::HTTPMethod::POST)
    ([](const crow::request& req) {
        crow::json::rvalue body;
        try
        {
            body = crow::json::load(req.body);
            if (!body)
            {
                return crow::response(400, R"({"error":"Invalid JSON"})");
            }
        }
        catch (...)
        {
            return crow::response(400, R"({"error":"Could not parse JSON"})");
        }

        // ==========================
        // 1. Parse the input JSON
        // ==========================
        auto planner_input = body["planner_input"];
        if (!planner_input)
        {
            return crow::response(400, R"({"error":"Missing 'planner_input' JSON block."})");
        }

        std::string map_file;
        try
        {
            map_file = planner_input["map_file"].s();
        }
        catch (...)
        {
            return crow::response(400, R"({"error":"Missing or invalid map_file property."})");
        }

        // Extract start pose
        double start_x, start_y, start_z, start_yaw;
        try
        {
            start_x   = planner_input["start_pose"]["x"].d();
            start_y   = planner_input["start_pose"]["y"].d();
            start_z   = planner_input["start_pose"]["z"].d();
            start_yaw = planner_input["start_pose"]["yaw"].d();
        }
        catch (...)
        {
            return crow::response(400, R"({"error":"Missing or invalid start_pose {x,y,z,yaw}."})");
        }

        // Extract goal region
        double goal_x, goal_y, goal_z, goal_radius;
        try
        {
            goal_x      = planner_input["goal_region"]["x"].d();
            goal_y      = planner_input["goal_region"]["y"].d();
            goal_z      = planner_input["goal_region"]["z"].d();
            goal_radius = planner_input["goal_region"]["radius"].d();
        }
        catch (...)
        {
            return crow::response(400, R"({"error":"Missing or invalid goal_region {x,y,z,radius}."})");
        }

        // Extract planner parameters
        double time_budget = 5.0;
        bool check_collision_max_altitude = true;
        double min_turn_radius = 60.0;
        double max_climb_angle = 0.15;
        try
        {
            time_budget = planner_input["planner_parameters"]["time_budget_seconds"].d();
            check_collision_max_altitude = planner_input["planner_parameters"]["check_collision_max_altitude"].b();
            min_turn_radius = planner_input["planner_parameters"]["dubins_min_turn_radius"].d();
            max_climb_angle = planner_input["planner_parameters"]["dubins_max_climb_angle"].d();
        }
        catch (...)
        {
            // Use defaults if missing
        }

        // ==========================
        // 2. Load map and setup planner
        // ==========================
        auto map_ptr = std::make_shared<TerrainMap>();
        bool success = false;
        try
        {
            success = map_ptr->initializeFromGeotiff(map_file);
            if (!success)
            {
                return crow::response(400, R"({"error":"Could not load the geotiff map."})");
            }
        }
        catch (const std::exception& e)
        {
            std::ostringstream oss;
            oss << R"({"error":"Exception while loading map: )" << e.what() << "\"}";
            return crow::response(400, oss.str());
        }

        // Create and configure the planner
        TerrainOmplRrt planner;
        planner.setMap(map_ptr);
        planner.setMaxAltitudeCollisionChecks(check_collision_max_altitude);
        planner.setAltitudeLimits(120.0, 50.0);  // Adjust limits as needed
        planner.setBoundsFromMap(map_ptr->getGridMap());

        planner.getProblemSetup()->getGeometricComponentStateSpace()
            ->as<fw_planning::spaces::DubinsAirplaneStateSpace>()->setMinTurningRadius(min_turn_radius);
        planner.getProblemSetup()->getGeometricComponentStateSpace()
            ->as<fw_planning::spaces::DubinsAirplaneStateSpace>()->setMaxClimbingAngle(max_climb_angle);

        // Configure start state with derived velocity
        double vx = std::cos(start_yaw) * 10.0;
        double vy = std::sin(start_yaw) * 10.0;
        double vz = 0.0;
        Eigen::Vector3d start_pos(start_x, start_y, start_z);
        Eigen::Vector3d start_vel(vx, vy, vz);

        // Setup goal region
        Eigen::Vector3d goal_pos(goal_x, goal_y, goal_z);
        planner.setupProblem(start_pos, start_vel, goal_pos, goal_radius);

        // ==========================
        // 3. Solve
        // ==========================
        std::vector<Eigen::Vector3d> solutionPath;
        bool found = planner.Solve(time_budget, solutionPath);
        double solve_time = planner.getSolutionTime();

        // Compute yaws for each segment
        std::vector<double> solutionYaws;
        solutionYaws.reserve(solutionPath.size());
        for (size_t i = 0; i < solutionPath.size(); i++)
        {
            if (i + 1 < solutionPath.size())
            {
                double yaw = std::atan2(solutionPath[i+1].y() - solutionPath[i].y(),
                                        solutionPath[i+1].x() - solutionPath[i].x());
                solutionYaws.push_back(yaw);
            }
            else
            {
                if (solutionYaws.empty())
                    solutionYaws.push_back(start_yaw);
                else
                    solutionYaws.push_back(solutionYaws.back());
            }
        }

        // ==========================
        // 4. Prepare JSON output
        // ==========================
        crow::json::wvalue out;
        out["planner_output"]["found_solution"] = found;
        out["planner_output"]["computation_time_seconds"] = solve_time;

        if (found && !solutionPath.empty()) {
            // Use move assignment to avoid copying the move-only wvalue
            crow::json::wvalue path_json = pathToJSON(solutionPath, solutionYaws);
            out["planner_output"]["path"] = std::move(path_json);
        } else {
            out["planner_output"]["path"] = crow::json::wvalue::list();
        }

        return crow::response{ out };
    });

    // Start the server on port 18080
    const uint16_t port = 18080;
    std::cout << "[TerrainPlannerServer] Listening on port " << port << "\n";
    app.port(port).multithreaded().run();

    return 0;
}