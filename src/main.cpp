//
// Created by maart on 5-6-2026.
//
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <limits>
#include <cmath>
#include <cassert>
#include <queue>
#include <tuple>
#include <vector>
#include <map>
#include <filesystem>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Point_3.h>
#include <CGAL/Segment_3.h>
#include <CGAL/intersections.h>

#include "json.hpp"

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Segment_3<Kernel> Segment;
typedef CGAL::Point_3<Kernel> Point;

const std::string input_file = "../input/Wellness_center.obj";
const std::string output_file = "../output/Wellness_center.city.json";

struct Shell {
    std::vector<Kernel::Triangle_3> triangles;
};

struct Object {
    std::string id;
    std::vector<Shell> shells;
};

std::map<std::string, Object> objects;

struct VoxelGrid {
    std::vector<unsigned int> voxels;
    unsigned int max_x, max_y, max_z;

    VoxelGrid(unsigned int x, unsigned int y, unsigned int z) {
        max_x = x;
        max_y = y;
        max_z = z;
        unsigned int total_voxels = x*y*z;
        voxels.reserve(total_voxels);
        for (unsigned int i = 0; i < total_voxels; ++i) voxels.push_back(0);
    }

    unsigned int &operator()(const unsigned int &x, const unsigned int &y, const unsigned int &z) {
        assert(x < max_x);
        assert(y < max_y);
        assert(z < max_z);
        return voxels[x + y*max_x + z*max_x*max_y];
    }

    unsigned int operator()(const unsigned int &x, const unsigned int &y, const unsigned int &z) const {
        assert(x < max_x);
        assert(y < max_y);
        assert(z < max_z);
        return voxels[x + y*max_x + z*max_x*max_y];
    }
};

void world_to_voxel(double wx, double wy, double wz,
                    double origin_x, double origin_y, double origin_z, double voxel_size,
                    unsigned int &i, unsigned int &j, unsigned int &k) {
    int calculated_i = static_cast<int>(std::floor((wx - origin_x) / voxel_size));
    int calculated_j = static_cast<int>(std::floor((wy - origin_y) / voxel_size));
    int calculated_k = static_cast<int>(std::floor((wz - origin_z) / voxel_size));

    i = static_cast<unsigned int>(std::max(0, calculated_i));
    j = static_cast<unsigned int>(std::max(0, calculated_j));
    k = static_cast<unsigned int>(std::max(0, calculated_k));
}

void voxel_to_world(unsigned int i, unsigned int j, unsigned int k,
                    double origin_x, double origin_y, double origin_z, double voxel_size,
                    double &wx, double &wy, double &wz) {

    wx = origin_x + (static_cast<double>(i)) * voxel_size;
    wy = origin_y + (static_cast<double>(j)) * voxel_size;
    wz = origin_z + (static_cast<double>(k)) * voxel_size;
}

int mark_voxel(VoxelGrid &voxel_grid, std::map<std::string, Object> objects,
    double origin_x, double origin_y, double origin_z, double voxel_size) {
    int marked = 0;
    for (const auto& [id, object] : objects) {
        for (const auto& shell : object.shells) {
            for (const auto& triangle : shell.triangles) {
                CGAL::Bbox_3 bbox = triangle.bbox();

                unsigned int min_i, min_j, min_k;
                unsigned int max_i, max_j, max_k;

                world_to_voxel(bbox.xmin(), bbox.ymin(), bbox.zmin(),
                    origin_x, origin_y, origin_z, voxel_size,
                    min_i, min_j, min_k
                    );

                world_to_voxel(bbox.xmax(), bbox.ymax(), bbox.zmax(),
                     origin_x, origin_y, origin_z, voxel_size,
                     max_i, max_j, max_k
                     );

                // Zorg dat de bounding box binnen de gridgrenzen blijft
                max_i = std::min(max_i, voxel_grid.max_x - 1);
                max_j = std::min(max_j, voxel_grid.max_y - 1);
                max_k = std::min(max_k, voxel_grid.max_z - 1);

                for (unsigned int i = min_i; i <= max_i; ++i) {
                    for (unsigned int j = min_j; j <= max_j; ++j) {
                        for (unsigned int k = min_k; k <= max_k; ++k) {

                            // Performance optimalisatie: check eerst of voxel al gevuld is
                            if (voxel_grid(i,j,k) == 1) continue;

                            double wx_min, wy_min, wz_min;
                            voxel_to_world(i, j, k,
                                origin_x, origin_y, origin_z, voxel_size,
                                wx_min, wy_min, wz_min
                                );

                            double wx_max = wx_min + voxel_size;
                            double wy_max = wy_min + voxel_size;
                            double wz_max = wz_min + voxel_size;

                            Point p1(wx_min, wy_min, wz_min);
                            Point p2(wx_max, wy_min, wz_min);
                            Point p3(wx_min, wy_max, wz_min);
                            Point p4(wx_max, wy_max, wz_min);
                            Point p5(wx_min, wy_min, wz_max);
                            Point p6(wx_max, wy_min, wz_max);
                            Point p7(wx_min, wy_max, wz_max);
                            Point p8(wx_max, wy_max, wz_max);

                            std::vector<Segment> edges;
                            edges.emplace_back(p1,p2); edges.emplace_back(p1,p3); edges.emplace_back(p1,p5);
                            edges.emplace_back(p2,p6); edges.emplace_back(p2,p4); edges.emplace_back(p3,p4);
                            edges.emplace_back(p3,p7); edges.emplace_back(p4,p8); edges.emplace_back(p5,p6);
                            edges.emplace_back(p5,p7); edges.emplace_back(p6,p8); edges.emplace_back(p7,p8);

                            for (const auto& edge: edges) {
                                if (CGAL::do_intersect(triangle, edge)) {
                                    voxel_grid(i,j,k) = 1;
                                    marked++;
                                    break;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return marked;
}

void mark_segment(VoxelGrid& voxel_grid, std::tuple<unsigned int, unsigned int, unsigned int> point, int segmentation_number, std::map<unsigned int, unsigned int>& segment_voxel_counter) {
    std::queue<std::tuple<unsigned int, unsigned int, unsigned int>> queue;

    unsigned int start_i = std::get<0>(point);
    unsigned int start_j = std::get<1>(point);
    unsigned int start_k = std::get<2>(point);

    if (voxel_grid(start_i, start_j, start_k) == 0) {
        voxel_grid(start_i, start_j, start_k) = segmentation_number;
        segment_voxel_counter[segmentation_number]++;
        queue.push(point);
    }

    const int dx[6] = { 1, -1, 0, 0, 0, 0 };
    const int dy[6] = { 0, 0, 1, -1, 0, 0 };
    const int dz[6] = { 0, 0, 0, 0, 1, -1 };

    while (!queue.empty()) {
        auto [i, j, k] = queue.front();
        queue.pop();

        for (int n = 0; n < 6; ++n) {
            int ni = static_cast<int>(i) + dx[n];
            int nj = static_cast<int>(j) + dy[n];
            int nk = static_cast<int>(k) + dz[n];

            if (ni < 0 || nj < 0 || nk < 0) continue;
            if (ni >= static_cast<int>(voxel_grid.max_x) ||
                nj >= static_cast<int>(voxel_grid.max_y) ||
                nk >= static_cast<int>(voxel_grid.max_z)) continue;

            if (voxel_grid(ni, nj, nk) == 0) {
                voxel_grid(ni, nj, nk) = segmentation_number;
                segment_voxel_counter[segmentation_number]++;
                queue.push({
                    static_cast<unsigned int>(ni),
                    static_cast<unsigned int>(nj),
                    static_cast<unsigned int>(nk)
                });
            }
        }
    }
}

void extract_surface(VoxelGrid& voxel_grid, double origin_x, double origin_y, double origin_z, double voxel_size,
    std::vector<std::vector<size_t>> &building_envelope_faces, std::map<unsigned int, std::vector<std::vector<size_t>>> &room_solid_faces, std::vector<std::array<double, 3>> &global_vertices) {

    std::cout << "\nStarten met Stap 6: Oppervlakken extraheren..." << std::endl;

    std::map<std::tuple<unsigned int, unsigned int, unsigned int>, size_t> corner_to_idx;

    // Function to save the corner coordinates
    auto get_vertex = [&](unsigned int u, unsigned int v, unsigned int w) -> size_t {
        std::tuple<unsigned int, unsigned int, unsigned int> corner = {u, v, w};
        if (corner_to_idx.count(corner)) return corner_to_idx[corner];

        double wx, wy, wz;
        voxel_to_world(u, v, w, origin_x, origin_y, origin_z, voxel_size, wx, wy, wz);
        global_vertices.push_back({wx, wy, wz});

        size_t new_idx = global_vertices.size() - 1;
        corner_to_idx[corner] = new_idx;
        return new_idx;
    };

    // Function to get the correct face
    auto get_face_indices = [&](unsigned int i, unsigned int j, unsigned int k, int n) -> std::vector<size_t> {
        std::vector<size_t> face;
        if (n == 0) { // +X (Oost)
            face.push_back(get_vertex(i + 1, j,     k));
            face.push_back(get_vertex(i + 1, j + 1, k));
            face.push_back(get_vertex(i + 1, j + 1, k + 1));
            face.push_back(get_vertex(i + 1, j,     k + 1));
        } else if (n == 1) { // -X (West)
            face.push_back(get_vertex(i,     j,     k));
            face.push_back(get_vertex(i,     j,     k + 1));
            face.push_back(get_vertex(i,     j + 1, k + 1));
            face.push_back(get_vertex(i,     j + 1, k));
        } else if (n == 2) { // +Y (Noord)
            face.push_back(get_vertex(i,     j + 1, k));
            face.push_back(get_vertex(i,     j + 1, k + 1));
            face.push_back(get_vertex(i + 1, j + 1, k + 1));
            face.push_back(get_vertex(i + 1, j + 1, k));
        } else if (n == 3) { // -Y (Zuid)
            face.push_back(get_vertex(i,     j,     k));
            face.push_back(get_vertex(i + 1, j,     k));
            face.push_back(get_vertex(i + 1, j,     k + 1));
            face.push_back(get_vertex(i,     j,     k + 1));
        } else if (n == 4) { // +Z (Boven)
            face.push_back(get_vertex(i,     j,     k + 1));
            face.push_back(get_vertex(i + 1, j,     k + 1));
            face.push_back(get_vertex(i + 1, j + 1, k + 1));
            face.push_back(get_vertex(i,     j + 1, k + 1));
        } else if (n == 5) { // -Z (Onder)
            face.push_back(get_vertex(i,     j,     k));
            face.push_back(get_vertex(i,     j + 1, k));
            face.push_back(get_vertex(i + 1, j + 1, k));
            face.push_back(get_vertex(i + 1, j,     k));
        }
        return face;
    };



    const int dx[6] = { 1, -1, 0, 0, 0, 0 };
    const int dy[6] = { 0, 0, 1, -1, 0, 0 };
    const int dz[6] = { 0, 0, 0, 0, 1, -1 };

    for (unsigned int i = 0; i < voxel_grid.max_x; ++i) {
        for (unsigned int j = 0; j < voxel_grid.max_y; ++j) {
            for (unsigned int k = 0; k < voxel_grid.max_z; ++k) {

                unsigned int current_id = voxel_grid(i, j, k);

                if (current_id != 1) {
                    continue;
                }

                unsigned int neighbor_id;

                for (int n = 0; n < 6; ++n) {
                    int ni = static_cast<int>(i) + dx[n];
                    int nj = static_cast<int>(j) + dy[n];
                    int nk = static_cast<int>(k) + dz[n];

                    if (ni < 0 || nj < 0 || nk < 0 ||
                    ni >= static_cast<int>(voxel_grid.max_x) ||
                    nj >= static_cast<int>(voxel_grid.max_y) ||
                    nk >= static_cast<int>(voxel_grid.max_z)) {
                        neighbor_id = 2;

                    } else {
                        neighbor_id = voxel_grid(ni, nj, nk);
                    }

                    if (neighbor_id == 2) {
                        building_envelope_faces.push_back(get_face_indices(i, j, k, n));
                    }
                    else if (current_id != neighbor_id) {
                        room_solid_faces[current_id].push_back(get_face_indices(i, j, k, n));
                    }
                }
            }
        }
    }
}

int main() {
    std::ifstream input_stream(input_file);
    if (!input_stream.is_open()) {
        std::cerr << "Kan bestand niet openen: " << input_file << std::endl;
        return 1;
    }

    std::string line;
    std::vector<Kernel::Point_3> vertices;
    std::string current_object_id = "";
    Shell current_shell;

    double min_x = std::numeric_limits<double>::infinity();
    double min_y = std::numeric_limits<double>::infinity();
    double min_z = std::numeric_limits<double>::infinity();
    double max_x = -std::numeric_limits<double>::infinity();
    double max_y = -std::numeric_limits<double>::infinity();
    double max_z = -std::numeric_limits<double>::infinity();

    while (getline(input_stream, line)) {
        std::istringstream line_stream(line);
        std::string line_type;
        line_stream >> line_type;

        if (line_type == "g") {
            if (!current_object_id.empty() && !current_shell.triangles.empty()) {
                objects[current_object_id].id = current_object_id;
                objects[current_object_id].shells.push_back(current_shell);
            }
            current_shell = Shell();
            line_stream >> current_object_id;
        }
        if (line_type == "v") {
            double x, y, z;
            line_stream >> x >> y >> z;
            vertices.push_back(Kernel::Point_3(x, y, z));
            min_x = std::min(min_x, x); min_y = std::min(min_y, y); min_z = std::min(min_z, z);
            max_x = std::max(max_x, x); max_y = std::max(max_y, y); max_z = std::max(max_z, z);
        }
        if (line_type == "f") {
            int idx1, idx2, idx3;
            char slash;
            int vn1, vn2, vn3;
            line_stream >> idx1 >> slash >> slash >> vn1;
            line_stream >> idx2 >> slash >> slash >> vn2;
            line_stream >> idx3 >> slash >> slash >> vn3;

            Kernel::Point_3 p1 = vertices[idx1 - 1];
            Kernel::Point_3 p2 = vertices[idx2 - 1];
            Kernel::Point_3 p3 = vertices[idx3 - 1];
            current_shell.triangles.push_back(Kernel::Triangle_3(p1, p2, p3));
        }
    }
    if (!current_object_id.empty() && !current_shell.triangles.empty()) {
        objects[current_object_id].id = current_object_id;
        objects[current_object_id].shells.push_back(current_shell);
    }

    double voxel_size = 0.5;
    unsigned int rows_x = static_cast<unsigned int>(std::ceil((max_x - min_x) / voxel_size));
    unsigned int rows_y = static_cast<unsigned int>(std::ceil((max_y - min_y) / voxel_size));
    unsigned int rows_z = static_cast<unsigned int>(std::ceil((max_z - min_z) / voxel_size));

    // VoxelGrid aanmaken met padding (+2)
    VoxelGrid my_building_grid(rows_x + 2, rows_y + 2, rows_z + 2);
    std::cout << "Grid succesvol aangemaakt met grootte: " << rows_x + 2 << "x" << rows_y + 2 << "x" << rows_z + 2 << std::endl;

    // FIX: Verschuif de origin met 1 voxel_size om de padding correct te benutten rondom het model
    double origin_x = min_x - voxel_size;
    double origin_y = min_y - voxel_size;
    double origin_z = min_z - voxel_size;

    int marked = mark_voxel(my_building_grid, objects, origin_x, origin_y, origin_z, voxel_size);
    std::cout << "Aantal unieke voxels gemarkeerd met 1 (muren): " << marked << std::endl;

    std::map<unsigned int, unsigned int> segment_voxel_counter;
    int segmentation_number = 2;
    std::tuple<unsigned int, unsigned int, unsigned int> exterior_point = {0, 0, 0};
    mark_segment(my_building_grid, exterior_point, segmentation_number, segment_voxel_counter);

    for (unsigned int i = 0; i < my_building_grid.max_x; ++i) {
        for (unsigned int j = 0; j < my_building_grid.max_y; ++j) {
            for (unsigned int k = 0; k < my_building_grid.max_z; ++k) {
                if (my_building_grid(i,j,k) == 0) {
                    segmentation_number++;
                    std::tuple<unsigned int, unsigned int, unsigned int> interior_point = {i, j, k};
                    mark_segment(my_building_grid, interior_point, segmentation_number, segment_voxel_counter);
                }
            }
        }
    }

    for (const auto& [segment_id, voxel_count] : segment_voxel_counter) {
        std::cout << "Segment " << segment_id << ": " << voxel_count << " voxels\n";
    }

    std::vector<std::vector<size_t>> building_envelope_faces;
    std::map<unsigned int, std::vector<std::vector<size_t>>> room_solid_faces;
    std::vector<std::array<double, 3>> global_vertices;


    extract_surface(my_building_grid, origin_x, origin_y, origin_z, voxel_size, building_envelope_faces, room_solid_faces, global_vertices);

    std::cout << "Starting CityJSON..." << std::endl;

    nlohmann::json json;
    json["type"] = "CityJSON";
    json["version"] = "2.0";
    json["transform"] = nlohmann::json::object();
    json["transform"]["scale"] = nlohmann::json::array({voxel_size, voxel_size, voxel_size});
    json["transform"]["translate"] = nlohmann::json::array({origin_x, origin_y, origin_z});
    json["CityObjects"] = nlohmann::json::object();

    // 1. Voeg het hoofdgebouw toe (Building)
    nlohmann::json building_obj;
    building_obj["type"] = "Building";

    nlohmann::json building_geom;
    building_geom["type"] = "Solid";
    building_geom["lod"] = "1.0";
    nlohmann::json building_surfaces = nlohmann::json::array();
    for(const auto& face : building_envelope_faces) {
        building_surfaces.push_back(nlohmann::json::array({ face }));
    }
    building_geom["boundaries"] = nlohmann::json::array({ building_surfaces });
    building_obj["geometry"] = nlohmann::json::array({ building_geom });

    json["CityObjects"]["id_main_building"] = building_obj;

    // 2. Voeg alle kamers toe (BuildingRoom)
    for (const auto& [room_id, faces] : room_solid_faces) {
        nlohmann::json room_obj;
        room_obj["type"] = "BuildingRoom";
        room_obj["parents"] = nlohmann::json::array({ "id_main_building" }); // Koppel kamer aan gebouw

        nlohmann::json room_geom;
        room_geom["type"] = "Solid";
        room_geom["lod"] = "1.0";
        nlohmann::json room_surfaces = nlohmann::json::array();
        for(const auto& face : faces) {
            room_surfaces.push_back(nlohmann::json::array({ face }));
        }
        room_geom["boundaries"] = nlohmann::json::array({ room_surfaces });
        room_obj["geometry"] = nlohmann::json::array({ room_geom });

        std::string room_key = "Room_" + std::to_string(room_id);
        json["CityObjects"][room_key] = room_obj;
    }

    // 3. Converteer de global_vertices naar de JSON vertex-lijst
    nlohmann::json json_vertices = nlohmann::json::array();
    for (const auto& v : global_vertices) {
        json_vertices.push_back({v[0], v[1], v[2]});
    }
    json["vertices"] = json_vertices;

    // 4. Schrijf de boel daadwerkelijk weg naar disk
    std::filesystem::path p(output_file);
    if (p.has_parent_path()) {
        std::filesystem::create_directories(p.parent_path());
    }

    std::ofstream out_file(output_file);
    if (out_file.is_open()) {
        out_file << json.dump(2); // Indentatie van 2 spaties voor leesbaarheid
        out_file.close();
        std::cout << "Hoera! CityJSON succesvol gegenereerd op: " << output_file << std::endl;
    } else {
        std::cerr << "Fout: Kon het CityJSON-outputbestand niet openen!" << std::endl;
    }

    return 0;
}