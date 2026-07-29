////////////////////////////////////////////////////////////////////////////////
// C++ include
#include <iostream>
#include <string>
#include <vector>
#include <limits>
#include <fstream>
#include <algorithm>
#include <numeric>

// Utilities for the Assignment
#include "utils.h"

// Image writing library
#define STB_IMAGE_WRITE_IMPLEMENTATION // Do not include this line twice in your project!
#include "stb_image_write.h"

// Shortcut to avoid Eigen:: everywhere, DO NOT USE IN .h
using namespace Eigen;

////////////////////////////////////////////////////////////////////////////////
// Class to store tree
////////////////////////////////////////////////////////////////////////////////
class AABBTree
{
public:
    class Node
    {
    public:
        AlignedBox3d bbox;
        int parent;   // Index of the parent node (-1 for root)
        int left;     // Index of the left child (-1 for a leaf)
        int right;    // Index of the right child (-1 for a leaf)
        int triangle; // Index of the node triangle (-1 for internal nodes)
    };

    std::vector<Node> nodes;
    int root;

    AABBTree() = default;                           // Default empty constructor
    AABBTree(const MatrixXd &V, const MatrixXi &F); // Build a BVH from an existing mesh

private:
    // builds the bvh recursively
    int build_recursive(const MatrixXd &V, const MatrixXi &F, const MatrixXd &centroids, int from, int to, int parent, std::vector<int> &triangles);
};

////////////////////////////////////////////////////////////////////////////////
// Scene setup, global variables
////////////////////////////////////////////////////////////////////////////////
const std::string data_dir = DATA_DIR;
const std::string filename("raytrace.png");
const std::string mesh_filename(data_dir + "bunny.off");

//Camera settings
const double focal_length = 2;
const double field_of_view = 0.7854; //45 degrees
const bool is_perspective = true;
const Vector3d camera_position(0, 0, 2);

// Triangle Mesh
MatrixXd vertices; // n x 3 matrix (n points)
MatrixXi facets;   // m x 3 matrix (m triangles)
AABBTree bvh;

//Material for the object, same material for all objects
const Vector4d obj_ambient_color(0.0, 0.5, 0.0, 0);
const Vector4d obj_diffuse_color(0.5, 0.5, 0.5, 0);
const Vector4d obj_specular_color(0.2, 0.2, 0.2, 0);
const double obj_specular_exponent = 256.0;
const Vector4d obj_reflection_color(0.7, 0.7, 0.7, 0);

// Precomputed (or otherwise) gradient vectors at each grid node
const int grid_size = 20;
std::vector<std::vector<Vector2d>> grid;

//Lights
std::vector<Vector3d> light_positions;
std::vector<Vector4d> light_colors;
//Ambient light
const Vector4d ambient_light(0.2, 0.2, 0.2, 0);

//Fills the different arrays
void setup_scene()
{
    //Loads file
    std::ifstream in(mesh_filename);
    std::string token;
    in >> token;
    int nv, nf, ne;
    in >> nv >> nf >> ne;
    vertices.resize(nv, 3);
    facets.resize(nf, 3);
    for (int i = 0; i < nv; ++i)
    {
        in >> vertices(i, 0) >> vertices(i, 1) >> vertices(i, 2);
    }
    for (int i = 0; i < nf; ++i)
    {
        int s;
        in >> s >> facets(i, 0) >> facets(i, 1) >> facets(i, 2);
        assert(s == 3);
    }

    //setup tree
    bvh = AABBTree(vertices, facets);

    //Lights
    light_positions.emplace_back(8, 8, 0);
    light_colors.emplace_back(16, 16, 16, 0);

    light_positions.emplace_back(6, -8, 0);
    light_colors.emplace_back(16, 16, 16, 0);

    light_positions.emplace_back(4, 8, 0);
    light_colors.emplace_back(16, 16, 16, 0);

    light_positions.emplace_back(2, -8, 0);
    light_colors.emplace_back(16, 16, 16, 0);

    light_positions.emplace_back(0, 8, 0);
    light_colors.emplace_back(16, 16, 16, 0);

    light_positions.emplace_back(-2, -8, 0);
    light_colors.emplace_back(16, 16, 16, 0);

    light_positions.emplace_back(-4, 8, 0);
    light_colors.emplace_back(16, 16, 16, 0);
}

////////////////////////////////////////////////////////////////////////////////
// BVH Code
////////////////////////////////////////////////////////////////////////////////

AlignedBox3d bbox_from_triangle(const Vector3d &a, const Vector3d &b, const Vector3d &c)
{
    AlignedBox3d box;
    box.extend(a);
    box.extend(b);
    box.extend(c);
    return box;
}

AABBTree::AABBTree(const MatrixXd &V, const MatrixXi &F)
{
    // Compute the centroids of all the triangles in the input mesh
    MatrixXd centroids(F.rows(), V.cols());
    centroids.setZero();
    for (int i = 0; i < F.rows(); ++i)
    {
        for (int k = 0; k < F.cols(); ++k)
        {
            centroids.row(i) += V.row(F(i, k));
        }
        centroids.row(i) /= F.cols();
    }

    //Vector containing the list of tringle indices
    std::vector<int> triangles(F.rows());
    std::iota(triangles.begin(), triangles.end(), 0);

    root = build_recursive(V, F, centroids, 0, triangles.size(), -1, triangles);
}

int AABBTree::build_recursive(const MatrixXd &V, const MatrixXi &F, const MatrixXd &centroids, int from, int to, int parent, std::vector<int> &triangles)
{
    // Scene is empty, so is the aabb tree
    if (to - from == 0)
    {
        return -1;
    }

    // If there is only 1 triangle left, then we are at a leaf
    if (to - from == 1)
    {
        int tri = triangles[from];
        int leaf = nodes.size();
        nodes.push_back(AABBTree::Node());
        Vector3d a = V.row(F(tri, 0));
        Vector3d b = V.row(F(tri, 1));
        Vector3d c = V.row(F(tri, 2));
        nodes[leaf].bbox = bbox_from_triangle(a,b,c);
        nodes[leaf].parent = parent;
        nodes[leaf].left = -1;
        nodes[leaf].right = -1;
        nodes[leaf].triangle = tri;
        return leaf;
    }

    // TODO sort centroids along the longest dimension

    //TODO Use AlignedBox3d to find the box around the current centroids
    AlignedBox3d centroid_box;
    centroid_box.setEmpty();
    for (int i = from; i < to; ++i) {
        int tri = triangles[i];
        Vector3d c = centroids.row(tri).transpose();
        centroid_box.extend(c);
    }

    Vector3d extent = centroid_box.diagonal();
    int longest_dim = 0;
    if (extent[1] > extent[longest_dim]) 
        longest_dim = 1;
    if (extent[2] > extent[longest_dim]) 
        longest_dim = 2;

    std::sort(triangles.begin() + from, triangles.begin() + to, [&](int f1, int f2) {
        return centroids(f1, longest_dim) < centroids(f2, longest_dim);
    });

    //TODO Create a n internal node and do a recursive call to build the left and right part of the tree
    int mid = from + (to - from) / 2;
    int node_idx = nodes.size();
    nodes.push_back(AABBTree::Node());
    nodes[node_idx].parent = parent;
    nodes[node_idx].triangle = -1;

    int left = build_recursive(
        V,
        F,
        centroids,
        from,
        mid,
        node_idx,
        triangles
    );
    int right = build_recursive(
        V,
        F,
        centroids,
        mid,
        to,
        node_idx,
        triangles
    );

    nodes[node_idx].left = left;
    nodes[node_idx].right = right;

    nodes[node_idx].bbox.setEmpty();
    if (left != -1) {
        nodes[node_idx].bbox.extend(nodes[left].bbox.min());
        nodes[node_idx].bbox.extend(nodes[left].bbox.max());
    }

    if (right != -1) {
        nodes[node_idx].bbox.extend(nodes[right].bbox.min());
        nodes[node_idx].bbox.extend(nodes[right].bbox.max());
    }

    return node_idx;
}

////////////////////////////////////////////////////////////////////////////////
// Intersection code
////////////////////////////////////////////////////////////////////////////////

double ray_triangle_intersection(const Vector3d &ray_origin, const Vector3d &ray_direction, const Vector3d &a, const Vector3d &b, const Vector3d &c, Vector3d &p, Vector3d &N)
{
    Vector3d e1 = b - a;
    Vector3d e2 = c - a;

    Matrix3d A;

    A.col(0) = e1;
    A.col(1) = e2;
    A.col(2) = -ray_direction;

    if (std::abs(A.determinant()) < 1e-10)
        return -1;

    Vector3d rhs = ray_origin - a;

    Vector3d x = A.colPivHouseholderQr().solve(rhs);

    double u = x(0);
    double v = x(1);
    double t = x(2);

    if (t < 0)
        return -1;

    if (u < 0 || v < 0 || u + v > 1)
        return -1;
    p = ray_origin + t * ray_direction;

    N = e1.cross(e2).normalized();

    return t;
}

bool ray_box_intersection(const Vector3d &ray_origin, const Vector3d &ray_direction, const AlignedBox3d &box)
{
    // TODO
    // Compute whether the ray intersects the given box.
    // we are not testing with the real surface here anyway.
    double tmin = -std::numeric_limits<double>::infinity();
    double tmax =  std::numeric_limits<double>::infinity();

    Vector3d bmin = box.min();
    Vector3d bmax = box.max();

    for (int i = 0; i < 3; i++)
    {
        if (std::abs(ray_direction[i]) < 1e-12)
        {
            if (ray_origin[i] < bmin[i] || ray_origin[i] > bmax[i])
                return false;
        }
        else
        {
            double t1 = (bmin[i] - ray_origin[i]) / ray_direction[i];
            double t2 = (bmax[i] - ray_origin[i]) / ray_direction[i];

            if (t1 > t2)
                std::swap(t1, t2);

            tmin = std::max(tmin, t1);
            tmax = std::min(tmax, t2);

            if (tmin > tmax)
                return false;
        }
    }

    return tmax >= std::max(tmin, 0.0);
}

//Finds the closest intersecting object returns its index
//In case of intersection it writes into p and N (intersection point and normals)
bool find_nearest_object(const Vector3d &ray_origin, const Vector3d &ray_direction, Vector3d &p, Vector3d &N)
{
    Vector3d tmp_p, tmp_N;
    double closest_t = std::numeric_limits<double>::max();
    bool intersects = false;
    // TODO
    // Method (1): Traverse every triangle and return the closest hit.
    // Method (2): Traverse the BVH tree and test the intersection with a
    // triangles at the leaf nodes that intersects the input ray.
    // for (int t = 0; t < facets.rows(); ++t) {
    //     int i0 = facets(t, 0);
    //     int i1 = facets(t, 1);
    //     int i2 = facets(t, 2);

    //     Vector3d a = vertices.row(i0);
    //     Vector3d b = vertices.row(i1);
    //     Vector3d c = vertices.row(i2);

    //     double result = ray_triangle_intersection(
    //         ray_origin,
    //         ray_direction,
    //         a,
    //         b,
    //         c,
    //         tmp_p,
    //         tmp_N
    //     );
    //     if (result != -1) {
    //         intersects = true;
    //         if (result < closest_t)
    //             closest_t = result;
    //     }
    //     p = tmp_p;
    //     N = tmp_N;
    // }
    std::vector<int> stack;

    if (bvh.root != -1)
        stack.push_back(bvh.root);

    while (!stack.empty())
    {
        int node_idx = stack.back();
        stack.pop_back();

        const AABBTree::Node &node = bvh.nodes[node_idx];

        if (!ray_box_intersection(ray_origin, ray_direction, node.bbox))
            continue;

        // Leaf node
        if (node.left == -1 && node.right == -1)
        {
            int t = node.triangle;

            int i0 = facets(t, 0);
            int i1 = facets(t, 1);
            int i2 = facets(t, 2);

            Vector3d a = vertices.row(i0);
            Vector3d b = vertices.row(i1);
            Vector3d c = vertices.row(i2);

            double hit = ray_triangle_intersection(
                ray_origin,
                ray_direction,
                a,
                b,
                c,
                tmp_p,
                tmp_N);

            if (hit != -1 && hit < closest_t)
            {
                closest_t = hit;
                p = tmp_p;
                N = tmp_N;
                intersects = true;
            }
        }
        else
        {
            if (node.left != -1)
                stack.push_back(node.left);

            if (node.right != -1)
                stack.push_back(node.right);
        }
    }

    return intersects;
}

////////////////////////////////////////////////////////////////////////////////
// Raytracer code
////////////////////////////////////////////////////////////////////////////////

Vector4d shoot_ray(const Vector3d &ray_origin, const Vector3d &ray_direction)
{
    //Intersection point and normal, these are output of find_nearest_object
    Vector3d p, N;

    const bool nearest_object = find_nearest_object(ray_origin, ray_direction, p, N);

    if (!nearest_object)
    {
        // Return a transparent color
        return Vector4d(0, 0, 0, 0);
    }

    // Ensure normal faces the incoming ray (toward the camera)
    if (N.dot(-ray_direction) < 0)
        N = -N;

    // Ambient light contribution
    const Vector4d ambient_color = obj_ambient_color.array() * ambient_light.array();

    // Direct lighting with shadows
    Vector4d lights_color(0, 0, 0, 0);
    for (int i = 0; i < (int)light_positions.size(); ++i)
    {
        const Vector3d &light_position = light_positions[i];
        const Vector4d &light_color = light_colors[i];

        Vector3d to_light = light_position - p;
        double dist2 = to_light.squaredNorm();
        Vector3d Li = to_light.normalized();

        // Shadow ray: offset origin to avoid self-intersection
        Vector3d shadow_origin = p + 1e-6 * N;
        Vector3d shadow_p, shadow_N;
        bool occluded = find_nearest_object(shadow_origin, Li, shadow_p, shadow_N);
        if (occluded && (shadow_p - p).squaredNorm() < dist2 - 1e-9)
            continue; // light is occluded

        // Diffuse
        double diff = std::max(Li.dot(N), 0.0);
        Vector4d diffuse = obj_diffuse_color * diff;

        // Specular (Blinn-Phong half-vector)
        Vector3d view_dir = (-ray_direction).normalized();
        Vector3d H = (Li + view_dir).normalized();
        double spec_k = std::pow(std::max(N.dot(H), 0.0), obj_specular_exponent);
        Vector4d specular = obj_specular_color * spec_k;

        lights_color += (diffuse + specular).cwiseProduct(light_color) / dist2;
    }

    // Rendering equation
    Vector4d C = ambient_color + lights_color;
    C(3) = 1; // opaque
    return C;
}

////////////////////////////////////////////////////////////////////////////////

void raytrace_scene()
{
    std::cout << "Simple ray tracer." << std::endl;

    int w = 640;
    int h = 480;
    MatrixXd R = MatrixXd::Zero(w, h);
    MatrixXd G = MatrixXd::Zero(w, h);
    MatrixXd B = MatrixXd::Zero(w, h);
    MatrixXd A = MatrixXd::Zero(w, h); // Store the alpha mask

    // The camera always points in the direction -z
    // The sensor grid is at a distance 'focal_length' from the camera center,
    // and covers an viewing angle given by 'field_of_view'.
    double aspect_ratio = double(w) / double(h);
    //TODO
    double image_y = std::tan(field_of_view / 2.0) * focal_length;
    double image_x = image_y * aspect_ratio;

    // The pixel grid through which we shoot rays is at a distance 'focal_length'
    const Vector3d image_origin(-image_x, image_y, camera_position[2] - focal_length);
    const Vector3d x_displacement(2.0 / w * image_x, 0, 0);
    const Vector3d y_displacement(0, -2.0 / h * image_y, 0);

    for (unsigned i = 0; i < w; ++i)
    {
        for (unsigned j = 0; j < h; ++j)
        {
            const Vector3d pixel_center = image_origin + (i + 0.5) * x_displacement + (j + 0.5) * y_displacement;

            // Prepare the ray
            Vector3d ray_origin;
            Vector3d ray_direction;

            if (is_perspective)
            {
                // Perspective camera
                ray_origin = camera_position;
                ray_direction = (pixel_center - camera_position).normalized();
            }
            else
            {
                // Orthographic camera
                ray_origin = pixel_center;
                ray_direction = Vector3d(0, 0, -1);
            }

            const Vector4d C = shoot_ray(ray_origin, ray_direction);
            R(i, j) = C(0);
            G(i, j) = C(1);
            B(i, j) = C(2);
            A(i, j) = C(3);
        }
    }

    // Save to png
    write_matrix_to_png(R, G, B, A, filename);
}

////////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{
    setup_scene();

    raytrace_scene();
    return 0;
}
