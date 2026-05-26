////////////////////////////////////////////////////////////////////////////////
#include <algorithm>
#include <complex>
#include <fstream>
#include <iostream>
#include <numeric>
#include <vector>
#include <string>

#include <Eigen/Dense>
// Shortcut to avoid  everywhere, DO NOT USE IN .h
using namespace Eigen;
////////////////////////////////////////////////////////////////////////////////

const std::string root_path = DATA_DIR;

// Computes the determinant of the matrix whose columns are the vector u and v
double inline det(const Vector2d &u, const Vector2d &v)
{
    return u.x() * v.y() - u.y() * v.x();
}

// Return true iff [a,b] intersects [c,d]
bool intersect_segment(const Vector2d &a, const Vector2d &b, const Vector2d &c, const Vector2d &d)
{
    double d1 = det(b - a, c - a);
    double d2 = det(b - a, d - a);
    double d3 = det(d - c, a - c);
    double d4 = det(d - c, b - c);

    return (d1 * d2 <= 0) && (d3 * d4 <= 0);
}

////////////////////////////////////////////////////////////////////////////////

bool is_inside(const std::vector<Vector2d> &poly, const Vector2d &query)
{
    double max_x, max_y;
    max_x = max_y = std::numeric_limits<double>::lowest();
    for (const Vector2d &p : poly)    {
        max_x = std::max(max_x, p.x());
        max_y = std::max(max_y, p.y());
    }
    Vector2d far_point(max_x + 1, max_y + 1);
    int count = 0;
    for (size_t i = 0; i < poly.size(); ++i)
    {
        if (intersect_segment(poly[i], poly[(i + 1) % poly.size()], query, far_point))
        {
            count++;
        }
    }
    std::cout << "count: " << count << std::endl;
    if (count % 2 == 1)
    {
        return true;
    }
    else
    {
        return false;
    }
}

////////////////////////////////////////////////////////////////////////////////

std::vector<Vector2d> load_xyz(const std::string &filename)
{
    std::ifstream in(filename);
    std::string line;
    std::getline(in, line, '\n');
    int num_points = std::stoi(line);
    std::vector<Vector2d> points;
    double x, y, z; 
    while (in >> x >> y >> z)
    {
        points.push_back(Vector2d(x, y));
    }
    return points;
}

void save_xyz(const std::string &filename, const std::vector<Vector2d> &points)
{
    std::ofstream out(root_path + "/" + filename);
    out << points.size() << std::endl;
    for(const Vector2d& p : points)
    {
        out << p.x() << " " << p.y() << " 0" << std::endl;
    }
}

std::vector<Vector2d> load_obj(const std::string &filename)
{
    std::ifstream in(filename);
    std::vector<Vector2d> points;
    std::vector<Vector2d> poly;
    char key;
    while (in >> key)
    {
        if (key == 'v')
        {
            double x, y, z;
            in >> x >> y >> z;
            points.push_back(Vector2d(x, y));
        }
        else if (key == 'f')
        {
            std::string line;
            std::getline(in, line);
            std::istringstream ss(line);
            int id;
            while (ss >> id)
            {
                poly.push_back(points[id - 1]);
            }
        }
    }
    return poly;
}

////////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{
    const std::string points_path = root_path + "/points.xyz";
    const std::string poly_path = root_path + "/polygon.obj";

    std::vector<Vector2d> points = load_xyz(points_path);

    ////////////////////////////////////////////////////////////////////////////////
    //Point in polygon
    std::vector<Vector2d> poly = load_obj(poly_path);
    std::vector<Vector2d> result;
    for (size_t i = 0; i < points.size(); ++i)
    {
        if (is_inside(poly, points[i]))
        {
            result.push_back(points[i]);
        }
    }
    save_xyz("output.xyz", result);

    return 0;
}
