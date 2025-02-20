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


#include <cmath>


void GeoConversions_forward(const double lat, const double lon, const double alt, double &y, double &x, double &h) {
  // 1. Convert the ellipsoidal latitudes φ and longitudes λ into arcseconds ["]
  const double lat_arc = lat * 3600.0;
  const double lon_arc = lon * 3600.0;

  // 2. Calculate the auxiliary values (differences of latitude and longitude relative to Bern in the unit [10000"]):
  //  φ' = (φ – 169028.66 ")/10000
  //  λ' = (λ – 26782.5 ")/10000
  const double lat_aux = (lat_arc - 169028.66) / 10000.0;
  const double lon_aux = (lon_arc - 26782.5) / 10000.0;

  // 3. Calculate projection coordinates in LV95 (E, N, h) or in LV03 (y, x, h)
  // E [m] = 2600072.37 + 211455.93 * λ' - 10938.51 * λ' * φ' - 0.36 * λ' * φ'2 - 44.54 * λ'3
  // y [m] = E – 2000000.00 N [m] = 1200147.07 + 308807.95 * φ' + 3745.25 * λ'2 + 76.63 * φ'2 - 194.56 * λ'2 * φ' +
  // 119.79 * φ'3 x [m] = N – 1000000.00
  // hCH [m] =hWGS – 49.55 + 2.73 * λ' + 6.94 * φ'
  const double E = 2600072.37 + 211455.93 * lon_aux - 10938.51 * lon_aux * lat_aux -
                   0.36 * lon_aux * std::pow(lat_aux, 2) - 44.54 * std::pow(lon_aux, 3);
  y = E - 2000000.00;
  const double N = 1200147.07 + 308807.95 * lat_aux + 3745.25 * std::pow(lon_aux, 2) + 76.63 * std::pow(lat_aux, 2) -
                   194.56 * std::pow(lon_aux, 2) * lat_aux + 119.79 * std::pow(lat_aux, 3);
  x = N - 1000000.00;

  h = alt - 49.55 + 2.73 * lon_aux + 6.84 * lat_aux;
};

void GeoConversions_reverse(const double y, const double x, const double h, double &lat, double &lon, double &alt) {
  // 1. Convert the projection coordinates E (easting) and N (northing) in LV95 (or y / x in LV03) into the civilian
  // system (Bern = 0 / 0) and express in the unit [1000 km]: E' = (E – 2600000 m)/1000000 = (y – 600000 m)/1000000
  // N' = (N – 1200000 m)/1000000 = (x – 200000 m)/1000000
  const double y_aux = (y - 600000.0) / 1000000.0;
  const double x_aux = (x - 200000.0) / 1000000.0;

  // 2. Calculate longitude λ and latitude φ in the unit [10000"]:
  //  λ' = 2.6779094 + 4.728982 * y' + 0.791484* y' * x' + 0.1306 * y' * x'2 - 0.0436 * y'3
  //  φ' = 16.9023892 + 3.238272 * x' - 0.270978 * y'2 - 0.002528 * x'2 - 0.0447 * y'2 * x' - 0.0140 * x'3
  // hWGS [m] = hCH + 49.55 - 12.60 * y' - 22.64 * x'
  const double lon_aux = 2.6779094 + 4.728982 * y_aux + 0.791484 * y_aux * x_aux + 0.1306 * y_aux * std::pow(x_aux, 2) -
                         0.0436 * std::pow(y_aux, 3);
  const double lat_aux = 16.9023892 + 3.238272 * x_aux - 0.270978 * std::pow(y_aux, 2) - 0.002528 * std::pow(x_aux, 2) -
                         0.0447 * std::pow(y_aux, 2) * x_aux - 0.0140 * std::pow(x_aux, 3);
  alt = h + 49.55 - 12.60 * y_aux - 22.64 * x_aux;

  lon = lon_aux * 100.0 / 36.0;
  lat = lat_aux * 100.0 / 36.0;
};


// A helper to convert a planned path (list of 3D states + yaw) to a JSON array
static crow::json::wvalue pathToJSON(const std::vector<Eigen::Vector3d>& path_points,
                                     const std::vector<double>& path_yaws,
                                     const std::vector<double>& path_lats,
                                     const std::vector<double>& path_lons,
                                     const std::vector<double>& path_alts)
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
        pt["lat"] = (i < path_lats.size()) ? path_lats[i] : 0.0;
        pt["lon"] = (i < path_lons.size()) ? path_lons[i] : 0.0;
        pt["alt"] = (i < path_alts.size()) ? path_alts[i] : 0.0;
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

//void lobalOriginCallback(double lat, double lon, double alt) {
//    std::cout << "[TerrainPlanner] Received Global Origin from FMU" << std::endl;
//  
//    // receive geocentric LLA coordinate from mavros/global_position/gp_origin
//    double geocentric_lat = static_cast<double>(lat);
//    double geocentric_lon = static_cast<double>(lon);
//    double geocentric_alt = static_cast<double>(alt);
//  
//    // convert to geodetic coordinates with WGS-84 ellipsoid as datum
//    GeographicLib::Geocentric earth(GeographicLib::Constants::WGS84_a(), GeographicLib::Constants::WGS84_f());
//    double geodetic_lat, geodetic_lon, geodetic_alt;
//    earth.Reverse(geocentric_lat, geocentric_lon, geocentric_alt, geodetic_lat, geodetic_lon, geodetic_alt);
//  
//    // create local cartesian coordinates (ENU)
//    enu_.emplace(geodetic_lat, geodetic_lon, geodetic_alt, GeographicLib::Geocentric::WGS84());
//  
//    // store the geodetic coordinates (why is this called local origin?)
//    local_origin_altitude_ = geodetic_alt;
//    local_origin_latitude_ = geodetic_lat;
//    local_origin_longitude_ = geodetic_lon;
//  
//    std::cout << "Global Origin" << std::endl;
//    std::cout << "lat: " << local_origin_latitude_ << std::endl;
//    std::cout << "lon: " << local_origin_longitude_ << std::endl;
//    std::cout << "alt: " << local_origin_altitude_ << std::endl;
//}

//std::ostream& operator<<(std::ostream& os, const ESPG& coord) {
//    // Assuming ESPG has members x, y, z. Adjust as needed.
//    os << "(" << coord.x << ", " << coord.y << ", " << coord.z << ")";
//    return os;
//}

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
        double start_lat, start_lon, start_alt;
        double start_x, start_y, start_z, start_yaw;
        try
        {
            start_lat   = planner_input["start_pose"]["x"].d();
            start_lon   = planner_input["start_pose"]["y"].d();
            start_alt   = planner_input["start_pose"]["z"].d();
            start_yaw = planner_input["start_pose"]["yaw"].d();
        }
        catch (...)
        {
            return crow::response(400, R"({"error":"Missing or invalid start_pose {x,y,z,yaw}."})");
        }

        GeoConversions_forward(start_lat, start_lon, start_alt, start_x, start_y, start_z);
        std::cout << "Start position0: " << start_x << ", " << start_y << ", " << start_z << std::endl;

        // Extract goal region
        double goal_lat, goal_lon, goal_alt;
        double goal_x, goal_y, goal_z, goal_radius;
        try
        {
            goal_lat      = planner_input["goal_region"]["x"].d();
            goal_lon      = planner_input["goal_region"]["y"].d();
            goal_alt      = planner_input["goal_region"]["z"].d();
            goal_radius = planner_input["goal_region"]["radius"].d();
        }
        catch (...)
        {
            return crow::response(400, R"({"error":"Missing or invalid goal_region {x,y,z,radius}."})");
        }
        GeoConversions_forward(goal_lat, goal_lon, goal_alt, goal_x, goal_y, goal_z);
        std::cout << "Goal position: " << goal_x << ", " << goal_y << ", " << goal_z << std::endl;

        //GeoConversions_reverse(const double y, const double x, const double h, double &lat, double &lon, double &alt)

        // Extract planner parameters
        double time_budget = 5.0;
        bool check_collision_max_altitude = true;
        double min_turn_radius = 60.0;
        double max_climb_angle = 0.15;

        double local_origin_altitude = -6.35513e+06;
        double local_origin_latitude = 89.9384;
        double local_origin_longitude = 11.8799;
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
        Eigen::Vector3d map_origin;

        try
        {
            success = map_ptr->initializeFromGeotiff(map_file);
            if (!success)
            {
                return crow::response(400, R"({"error":"Could not load the geotiff map."})");
            }

            map_ptr->AddLayerDistanceTransform(50, "distance_surface");
            map_ptr->AddLayerDistanceTransform(120, "max_elevation");
            map_ptr->AddLayerHorizontalDistanceTransform(80, "ics_+", "distance_surface");
            map_ptr->AddLayerHorizontalDistanceTransform(-80, "ics_-", "max_elevation");
            map_ptr->addLayerSafety("safety", "ics_+", "ics_-");

            ESPG map_coordinate;
            map_ptr->getGlobalOrigin(map_coordinate, map_origin);
            start_x -= map_origin.x();
            start_y -= map_origin.y();
            start_z -= map_origin.z();

            std::cout << "map_origin: " << map_origin.x() << ", " << map_origin.y() << ", " << map_origin.z() << std::endl;
            std::cout << "Start position: " << start_x << ", " << start_y << ", " << start_z << std::endl;



            goal_x -= map_origin.x();
            goal_y -= map_origin.y();
            goal_z -= map_origin.z();
            std::cout << "Goal position: " << goal_x << ", " << goal_y << ", " << goal_z << std::endl;

            // For ESPG, ensure an operator<< is defined. Otherwise, print its members manually.
            //std::cout << "Map Coordinate: " << map_coordinate << std::endl;

            // Eigen vectors have an overloaded << operator.
            // By default, printing an Eigen::Vector3d will output each element on a separate line.
            // If you prefer a single-line output, use .transpose():
            std::cout << "Map Origin: " << map_origin.transpose() << std::endl;

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
        planner.setAltitudeLimits(4000.0, 0.0);
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
        // planner.setupProblem(start_pos, start_vel, goal_pos, goal_radius);
        std::cout << "Setting up problem with start: " << start_pos.transpose() << " and goal: " << goal_pos.transpose() <<  goal_radius << std::endl;
        planner.setupProblem(start_pos, goal_pos, goal_radius);

        // ==========================
        // 3. Solve
        // ==========================
        std::vector<Eigen::Vector3d> solutionPath;
        bool found = planner.Solve(time_budget, solutionPath);
        double solve_time = planner.getSolutionTime();

        // Compute yaws for each segment
        std::vector<double> solutionYaws;
        std::vector<double> solutionLats;
        std::vector<double> solutionLons;
        std::vector<double> solutionAlts;
        solutionYaws.reserve(solutionPath.size());
        for (size_t i = 0; i < solutionPath.size(); i++)
        {
            double lat, lon, alt;
            GeoConversions_reverse(solutionPath[i].x() + map_origin.x(), solutionPath[i].y() + map_origin.y(), solutionPath[i].z() + map_origin.z(), lat, lon, alt);
            solutionLats.push_back(lat);
            solutionLons.push_back(lon);
            solutionAlts.push_back(alt);

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
            crow::json::wvalue path_json = pathToJSON(solutionPath, solutionYaws, solutionLats, solutionLons, solutionAlts);
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