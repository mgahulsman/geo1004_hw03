//
// Created by maart on 5-6-2026.
//
#include <iostream>
#include <fstream>
#include <sstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <algorithm>
#include <limits>
#include <cmath>
#include <cassert>
#include <CGAL/Point_3.h>
#include <CGAL/Segment_3.h>
#include <CGAL/intersections.h>
#include <queue>
#include <tuple>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Segment_3<Kernel> Segment;
typedef CGAL::Point_3<Kernel> Point;

const std::string input_file = "../input/IfcOpenHouse_IFC4.obj";
const std::string output_file = "../output/IfcOpenHouse_IFC4.city.json";

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

    i = static_cast<unsigned int>(calculated_i);
    j = static_cast<unsigned int>(calculated_j);
    k = static_cast<unsigned int>(calculated_k);
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
    // loop through each traing
    int marked = 0;
    for (const auto& [id, object] : objects) {
        for (const auto& shell : object.shells) {
            for (const auto& triangle : shell.triangles) {
                // calc bbox for each triang
                CGAL::Bbox_3 bbox = triangle.bbox();

                // Get the min and max voxel coordinates of the triangle's bounding box
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

                // Loop through the voxel indices within the bounding box
                for (unsigned int i = min_i; i <= max_i; ++i) {
                    for (unsigned int j = min_j; j <= max_j; ++j) {
                        for (unsigned int k = min_k; k <= max_k; ++k) {
                            double wx_min, wy_min, wz_min;
                            voxel_to_world(i, j, k,
                                origin_x, origin_y, origin_z, voxel_size,
                                wx_min, wy_min, wz_min
                                );

                            double wx_max = wx_min + voxel_size;
                            double wy_max = wy_min + voxel_size;
                            double wz_max = wz_min + voxel_size;

                            // cgal do intersect ding
                            Point p1(wx_min, wy_min, wz_min);
                            Point p2(wx_max, wy_min, wz_min);
                            Point p3(wx_min, wy_max, wz_min);
                            Point p4(wx_max, wy_max, wz_min);
                            Point p5(wx_min, wy_min, wz_max);
                            Point p6(wx_max, wy_min, wz_max);
                            Point p7(wx_min, wy_max, wz_max);
                            Point p8(wx_max, wy_max, wz_max);

                            std::vector<Segment> edges;
                            edges.emplace_back(p1,p2);
                            edges.emplace_back(p1,p3);
                            edges.emplace_back(p1,p5);
                            edges.emplace_back(p2,p6);

                            edges.emplace_back(p2,p4);
                            edges.emplace_back(p3,p4);
                            edges.emplace_back(p3,p7);
                            edges.emplace_back(p4,p8);

                            edges.emplace_back(p5,p6);
                            edges.emplace_back(p5,p7);
                            edges.emplace_back(p6,p8);
                            edges.emplace_back(p7,p8);

                            for (const auto& edge: edges) {
                                if (CGAL::do_intersect(triangle, edge)) {
                                    voxel_grid(i,j,k) = 1;
                                    marked++;
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

void mark_exterior_voxel(VoxelGrid& voxel_grid)
{
    std::queue<std::tuple<unsigned int, unsigned int, unsigned int>> queue;

    queue.push({0, 0, 0});

    while (!queue.empty()) {
        auto [i, j, k] = queue.front();
        queue.pop();

        // Already visited?
        if (voxel_grid(i, j, k) != 0)
            continue;

        voxel_grid(i, j, k) = 2;

        const int dx[6] = { 1, -1, 0, 0, 0, 0 };
        const int dy[6] = { 0, 0, 1, -1, 0, 0 };
        const int dz[6] = { 0, 0, 0, 0, 1, -1 };

        for (int n = 0; n < 6; ++n)
        {
            int ni = static_cast<int>(i) + dx[n];
            int nj = static_cast<int>(j) + dy[n];
            int nk = static_cast<int>(k) + dz[n];

            // Bounds check
            if (ni < 0 || nj < 0 || nk < 0)
                continue;

            if (ni >= static_cast<int>(voxel_grid.max_x) ||
                nj >= static_cast<int>(voxel_grid.max_y) ||
                nk >= static_cast<int>(voxel_grid.max_z))
                continue;

            // Only flood-fill empty voxels
            if (voxel_grid(ni, nj, nk) == 0)
            {
                queue.push({
                    static_cast<unsigned int>(ni),
                    static_cast<unsigned int>(nj),
                    static_cast<unsigned int>(nk)
                });
            }
        }
    }
}

void mark_interior_voxel(VoxelGrid& voxel_grid, std::tuple<int, int, int> point, int segmentation_number) {
    std::queue<std::tuple<unsigned int, unsigned int, unsigned int>> queue;
    queue.push({point});

    while (!queue.empty()) {
        auto [i, j, k] = queue.front();
        queue.pop();

        // Already visited?
        if (voxel_grid(i, j, k) != 0)
            continue;

        voxel_grid(i, j, k) = segmentation_number;

        const int dx[6] = { 1, -1, 0, 0, 0, 0 };
        const int dy[6] = { 0, 0, 1, -1, 0, 0 };
        const int dz[6] = { 0, 0, 0, 0, 1, -1 };

        for (int n = 0; n < 6; ++n)
        {
            int ni = static_cast<int>(i) + dx[n];
            int nj = static_cast<int>(j) + dy[n];
            int nk = static_cast<int>(k) + dz[n];

            // // Bounds check - not necessary here because we are in an interior area
            // if (ni < 0 || nj < 0 || nk < 0)
            //     continue;
            //
            // if (ni >= static_cast<int>(voxel_grid.max_x) ||
            //     nj >= static_cast<int>(voxel_grid.max_y) ||
            //     nk >= static_cast<int>(voxel_grid.max_z))
            //     continue;

            // Only flood-fill empty voxels
            if (voxel_grid(ni, nj, nk) == 0)
            {
                queue.push({
                    static_cast<unsigned int>(ni),
                    static_cast<unsigned int>(nj),
                    static_cast<unsigned int>(nk)
                });
            }
        }
    }

}

int main() {
    std::ifstream input_stream;
    input_stream.open(input_file);

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

    // Parse line by line
    while (getline(input_stream, line)) {
        std::istringstream line_stream(line);
        std::string line_type;
        line_stream >> line_type;

        Shell shell;

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

            min_x = std::min(min_x, x);
            min_y = std::min(min_y, y);
            min_z = std::min(min_z, z);

            max_x = std::max(max_x, x);
            max_y = std::max(max_y, y);
            max_z = std::max(max_z, z);

        }

        if (line_type == "f") {
            int idx1, idx2, idx3;
            char slash;
            int vn1, vn2, vn3;

            // TODO: potential error: faces heeft ook een texture
            line_stream >> idx1 >> slash >> slash >> vn1;
            line_stream >> idx2 >> slash >> slash >> vn2;
            line_stream >> idx3 >> slash >> slash >> vn3;

            Kernel::Point_3 p1 = vertices[idx1 - 1];
            Kernel::Point_3 p2 = vertices[idx2 - 1];
            Kernel::Point_3 p3 = vertices[idx3 - 1];

            Kernel::Triangle_3 tri(p1, p2, p3);
            current_shell.triangles.push_back(tri);
        }
    }
    if (!current_object_id.empty() && !current_shell.triangles.empty()) {
        objects[current_object_id].id = current_object_id;
        objects[current_object_id].shells.push_back(current_shell);
    }
    else {
        std::cerr << "obj file empty?" << std::endl;
    }

    std::cout << "\nInlezen voltooid! Resultaten:" << std::endl;
    for (const auto& [id, obj] : objects) {
        std::cout << "Object: " << id << " heeft " << obj.shells.size() << " shell(s) met in totaal ";
        size_t total_triangles = 0;
        for (const auto& s : obj.shells) total_triangles += s.triangles.size();
        std::cout << total_triangles << " driehoeken." << std::endl;
    }

    double voxel_size = 0.5;
    unsigned int rows_x = static_cast<unsigned int>(std::ceil((max_x - min_x) / voxel_size));
    unsigned int rows_y = static_cast<unsigned int>(std::ceil((max_y - min_y) / voxel_size));
    unsigned int rows_z = static_cast<unsigned int>(std::ceil((max_z - min_z) / voxel_size));

    VoxelGrid my_building_grid(rows_x + 2, rows_y + 2, rows_z + 2);
    std::cout << "Grid succesvol aangemaakt met grootte: " << rows_x << "x" << rows_y << "x" << rows_z << std::endl;

    // Apply bbox heuristic
    int marked = mark_voxel(my_building_grid, objects, min_x, min_y, min_z, voxel_size);
    std::cout << "Aantal voxels gemarkeerd met 1: " << marked << std::endl;

    mark_exterior_voxel(my_building_grid);

    int interior_number = 3;
    for (unsigned int i = 0; i < my_building_grid.max_x; ++i)
    {
        for (unsigned int j = 0; j < my_building_grid.max_y; ++j)
        {
            for (unsigned int k = 0; k < my_building_grid.max_z; ++k)
            {
                if (my_building_grid(i,j,k) == 0)
                    std::tuple<int, int, int> interior_point = {i, j, k};
                    mark_interior_voxel(my_building_grid, interior_point, interior_number); // interior
                    interior_number++;
            }
        }
    }

    return 0;
}