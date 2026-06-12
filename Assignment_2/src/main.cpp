// C++ include
#include <iostream>
#include <string>
#include <vector>

// Utilities for the Assignment
#include "utils.h"

// Image writing library
#define STB_IMAGE_WRITE_IMPLEMENTATION // Do not include this line twice in your project!
#include "stb_image_write.h"

// Shortcut to avoid Eigen:: everywhere, DO NOT USE IN .h
using namespace Eigen;

const float EPSILON = 1e-6f;

void raytrace_sphere()
{
    std::cout << "Simple ray tracer, one sphere with orthographic projection" << std::endl;

    const std::string filename("sphere_orthographic.png");
    MatrixXd C = MatrixXd::Zero(800, 800); // Store the color
    MatrixXd A = MatrixXd::Zero(800, 800); // Store the alpha mask

    const Vector3d camera_origin(0, 0, 3);
    const Vector3d camera_view_direction(0, 0, -1);

    // The camera is orthographic, pointing in the direction -z and covering the
    // unit square (-1,1) in x and y
    const Vector3d image_origin(-1, 1, 1);
    const Vector3d x_displacement(2.0 / C.cols(), 0, 0);
    const Vector3d y_displacement(0, -2.0 / C.rows(), 0);

    const double sphere_radius = 0.9;
    const Vector3d sphere_center(0, 0, 0);

    // Single light source
    const Vector3d light_position(-1, 1, 1);

    for (unsigned i = 0; i < C.cols(); ++i)
    {
        for (unsigned j = 0; j < C.rows(); ++j)
        {
            const Vector3d pixel_center = image_origin + double(i) * x_displacement + double(j) * y_displacement;

            // Prepare the ray
            const Vector3d ray_origin = pixel_center;
            const Vector3d ray_direction = camera_view_direction;

            // Intersect with the sphere
            // NOTE: this is a special case of a sphere centered in the origin and for orthographic rays aligned with the z axis
            // TODO change this with the generic case
            Vector2d ray_on_xy(ray_origin(0), ray_origin(1));

            if (ray_on_xy.norm() < sphere_radius)
            {
                // The ray hit the sphere, compute the exact intersection point
                auto a = ray_direction.dot(ray_direction);
                auto b = 2 * ray_direction.dot(ray_origin - sphere_center);
                auto c = (ray_origin - sphere_center).dot(ray_origin - sphere_center) - sphere_radius * sphere_radius;
                auto discriminant = b * b - 4 * a * c;
                auto t1 = (-b - sqrt(discriminant)) / (2 * a);
                auto t2 = (-b + sqrt(discriminant)) / (2 * a);
                auto t = std::min(t1, t2);
                Vector3d ray_intersection = ray_origin + t * ray_direction;

                // Compute normal at the intersection point
                Vector3d ray_normal = ray_intersection.normalized();

                // Simple diffuse model
                C(i, j) = (light_position - ray_intersection).normalized().transpose() * ray_normal;

                // Clamp to zero
                C(i, j) = std::max(C(i, j), 0.);

                // Disable the alpha mask for this pixel
                A(i, j) = 1;
            }
        }
    }

    // Save to png
    write_matrix_to_png(C, C, C, A, filename);
}

void raytrace_parallelogram()
{
    std::cout << "Simple ray tracer, one parallelogram with orthographic projection" << std::endl;

    const std::string filename("plane_orthographic.png");
    MatrixXd C = MatrixXd::Zero(800, 800); // Store the color
    MatrixXd A = MatrixXd::Zero(800, 800); // Store the alpha mask

    const Vector3d camera_origin(0, 0, 3);
    const Vector3d camera_view_direction(0, 0, -1);

    // The camera is orthographic, pointing in the direction -z and covering the unit square (-1,1) in x and y
    const Vector3d image_origin(-1, 1, 1);
    const Vector3d x_displacement(2.0 / C.cols(), 0, 0);
    const Vector3d y_displacement(0, -2.0 / C.rows(), 0);

    // Parameters of the parallelogram (position of the lower-left corner + two sides)
    const Vector3d pgram_origin(-0.5, -0.5, 0);
    const Vector3d pgram_u(1, 0.4, 0);
    const Vector3d pgram_v(0, 0.7, -10);

    // Single light source
    const Vector3d light_position(-1, 1, 1);

    for (unsigned i = 0; i < C.cols(); ++i)
    {
        for (unsigned j = 0; j < C.rows(); ++j)
        {
            const Vector3d pixel_center = image_origin + double(i) * x_displacement + double(j) * y_displacement;

            // Prepare the ray
            const Vector3d ray_origin = pixel_center;
            const Vector3d ray_direction = camera_view_direction;

            // TODO: Check if the ray intersects with the parallelogram
            // TODO: The ray hit the parallelogram, compute the exact intersection
            // point
            auto w = ray_origin - pgram_origin;
            auto n = pgram_u.cross(pgram_v);

            if(std::fabs(n.dot(ray_direction)) < EPSILON)
                continue;

            auto distance_to_plane = n.dot(pgram_origin - ray_origin) / n.dot(ray_direction);

            if(distance_to_plane < 0)
                continue;

            auto ray_intersection = ray_origin + distance_to_plane * ray_direction;

            auto origin_to_intersection = ray_intersection - pgram_origin;

            auto s = origin_to_intersection.dot(pgram_u) * pgram_v.dot(pgram_v) - origin_to_intersection.dot(pgram_v) * pgram_u.dot(pgram_v);
            s /= pgram_u.dot(pgram_u) * pgram_v.dot(pgram_v) - pgram_u.dot(pgram_v) * pgram_u.dot(pgram_v);

            if(s < 0 || s > 1)
                continue;

            auto t = origin_to_intersection.dot(pgram_v) * pgram_u.dot(pgram_u) - origin_to_intersection.dot(pgram_u) * pgram_u.dot(pgram_v);
            t /= pgram_u.dot(pgram_u) * pgram_v.dot(pgram_v) - pgram_u.dot(pgram_v) * pgram_u.dot(pgram_v);

            if(t < 0 || t > 1)
                continue;


            // TODO: Compute normal at the intersection point
            auto ray_normal = n.normalized();

            // Simple diffuse model
            C(i, j) = (light_position - ray_intersection).normalized().transpose() * ray_normal;

            // Clamp to zero
            C(i, j) = std::max(C(i, j), 0.);

            // Disable the alpha mask for this pixel
            A(i, j) = 1;
        }
    }

    // Save to png
    write_matrix_to_png(C, C, C, A, filename);
}

void raytrace_perspective()
{
    std::cout << "Simple ray tracer, one parallelogram with perspective projection" << std::endl;

    const std::string filename("plane_perspective.png");
    MatrixXd C = MatrixXd::Zero(800, 800); // Store the color
    MatrixXd A = MatrixXd::Zero(800, 800); // Store the alpha mask

    const Vector3d camera_origin(0, 0, 3);
    const Vector3d camera_view_direction(0, 0, -1);

    // The camera is perspective, pointing in the direction -z and covering the unit square (-1,1) in x and y
    const Vector3d image_origin(-1, 1, 1);
    const Vector3d x_displacement(2.0 / C.cols(), 0, 0);
    const Vector3d y_displacement(0, -2.0 / C.rows(), 0);

    // TODO: Parameters of the parallelogram (position of the lower-left corner + two sides)
    const Vector3d pgram_origin(-0.5, -0.5, 0);
    const Vector3d pgram_u(1, 0.4, 0);
    const Vector3d pgram_v(0, 0.7, -10);

    // Single light source
    const Vector3d light_position(-1, 1, 1);

    for (unsigned i = 0; i < C.cols(); ++i)
    {
        for (unsigned j = 0; j < C.rows(); ++j)
        {
            const Vector3d pixel_center = image_origin + double(i) * x_displacement + double(j) * y_displacement;

            // TODO: Prepare the ray (origin point and direction)
            const Vector3d ray_origin = pixel_center;
            const Vector3d ray_direction = pixel_center - camera_origin;

            auto w = ray_origin - pgram_origin;
            auto n = pgram_u.cross(pgram_v);

            if(std::fabs(n.dot(ray_direction)) < EPSILON)
                continue;

            auto distance_to_plane = n.dot(pgram_origin - ray_origin) / n.dot(ray_direction);

            if(distance_to_plane < 0)
                continue;

            auto ray_intersection = ray_origin + distance_to_plane * ray_direction;

            auto origin_to_intersection = ray_intersection - pgram_origin;

            auto s = origin_to_intersection.dot(pgram_u) * pgram_v.dot(pgram_v) - origin_to_intersection.dot(pgram_v) * pgram_u.dot(pgram_v);
            s /= pgram_u.dot(pgram_u) * pgram_v.dot(pgram_v) - pgram_u.dot(pgram_v) * pgram_u.dot(pgram_v);

            if(s < 0 || s > 1)
                continue;

            auto t = origin_to_intersection.dot(pgram_v) * pgram_u.dot(pgram_u) - origin_to_intersection.dot(pgram_u) * pgram_u.dot(pgram_v);
            t /= pgram_u.dot(pgram_u) * pgram_v.dot(pgram_v) - pgram_u.dot(pgram_v) * pgram_u.dot(pgram_v);

            if(t < 0 || t > 1)
                continue;


            // TODO: Compute normal at the intersection point
            auto ray_normal = n.normalized();

            // Simple diffuse model
            C(i, j) = (light_position - ray_intersection).normalized().transpose() * ray_normal;

            // Clamp to zero
            C(i, j) = std::max(C(i, j), 0.);

            // Disable the alpha mask for this pixel
            A(i, j) = 1;
        }
    }

    // Save to png
    write_matrix_to_png(C, C, C, A, filename);
}

void raytrace_shading()
{
    std::cout << "Simple ray tracer, one sphere with different shading" << std::endl;

    const std::string filename("shading.png");
    Eigen::Matrix<Eigen::Vector3d, Eigen::Dynamic, Eigen::Dynamic> C(800, 800); // Store the color
    MatrixXd A = MatrixXd::Zero(800, 800); // Store the alpha mask

    const Vector3d camera_origin(0, 0, 3);
    const Vector3d camera_view_direction(0, 0, -1);

    // The camera is perspective, pointing in the direction -z and covering the unit square (-1,1) in x and y
    const Vector3d image_origin(-1, 1, 1);
    const Vector3d x_displacement(2.0 / A.cols(), 0, 0);
    const Vector3d y_displacement(0, -2.0 / A.rows(), 0);

    //Sphere setup
    const Vector3d sphere_center(0, 0, 0);
    const double sphere_radius = 0.9;

    //material params
    const Vector3d diffuse_color(1, 0, 1);
    const double specular_exponent = 100;
    const Vector3d specular_color(0., 0, 1);

    // Single light source
    const Vector3d light_position(-1, 1, 1);
    const Vector3d light_intesity(1, 1, 1);
    const Vector3d ambient = 0.1 * Vector3d(1, 1, 1);

    for (unsigned i = 0; i < C.cols(); ++i)
    {
        for (unsigned j = 0; j < C.rows(); ++j)
        {
            const Vector3d pixel_center = image_origin + double(i) * x_displacement + double(j) * y_displacement;

            // TODO: Prepare the ray (origin point and direction)
            const Vector3d ray_origin = pixel_center;
            const Vector3d ray_direction = camera_view_direction;

            Vector2d ray_on_xy(ray_origin(0), ray_origin(1));

            if (ray_on_xy.norm() < sphere_radius)
            {
                // The ray hit the sphere, compute the exact intersection point
                auto a = ray_direction.dot(ray_direction);
                auto b = 2 * ray_direction.dot(ray_origin - sphere_center);
                auto c = (ray_origin - sphere_center).dot(ray_origin - sphere_center) - sphere_radius * sphere_radius;
                auto discriminant = b * b - 4 * a * c;
                auto t1 = (-b - sqrt(discriminant)) / (2 * a);
                auto t2 = (-b + sqrt(discriminant)) / (2 * a);
                auto t = std::min(t1, t2);
                Vector3d ray_intersection = ray_origin + t * ray_direction;

                // Compute normal at the intersection point
                Vector3d ray_normal = ray_intersection.normalized();

                // TODO: Add shading parameter here
                auto L = (light_position - ray_intersection).normalized(); // light direction
                auto diffuse = diffuse_color * (std::max(L.dot(ray_normal), 0.0));
                auto V = -ray_direction.normalized(); // view direction
                auto R = (2 * ray_normal.dot(L) * ray_normal - L).normalized(); // reflection
                auto specular = specular_color * pow(std::max(R.dot(V), 0.0), specular_exponent);

                // Combine shading + ambient
                Vector3d pixel_color = ambient + diffuse + specular;

                // Clamp each channel to [0,1]
                pixel_color = pixel_color.cwiseMax(0.0).cwiseMin(1.0);

                // Assign to color matrix
                C(i, j) = pixel_color;

                // Set alpha
                A(i, j) = 1.0;
            }
        }
    }

    // Save to png
    write_matrix_to_png(C.unaryExpr([](const Vector3d &v) { return v(0); }),
        C.unaryExpr([](const Vector3d &v) { return v(1); }),
        C.unaryExpr([](const Vector3d &v) { return v(2); }), A, filename);
}

int main()
{
    raytrace_sphere();
    raytrace_parallelogram();
    raytrace_perspective();
    raytrace_shading();

    return 0;
}
