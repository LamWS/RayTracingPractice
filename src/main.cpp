#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include "Ray.h"
#include "Hitable.h"
#include <cfloat>
#include "Camera.h"
#include "random"
#include "Material.h"

using namespace std;
using namespace Eigen;


Vector3d color(const Ray &r, Hitable *world, int depth) {
    hit_record rec;
    if (world->hit(r, 0.001, FLT_MAX, rec)) {
        Ray scatter;
        Vector3d attenuation;
        if (depth < 50 && rec.material->scatter(r, rec, attenuation, scatter)) {
            auto c = color(scatter, world, depth + 1);
            return Vector3d(c.x() * attenuation.x(), c.y() * attenuation.x(), c.z() * attenuation.z());
        }
        return Vector3d(0, 0, 0);
    } else {
        Vector3d unit_direction = r.direction().normalized();
        double t = 0.5 * (unit_direction.y() + 1.0);
        return (1 - t) * Vector3d(1.0, 1.0, 1.0) + t * Vector3d(0.5, 0.7, 1.0);
    }
}

int main() {
    ofstream outputFile;
    outputFile.open("output.ppm", ios::out);

    int nx = 400, ny = 200, ns = 100;
    outputFile << "P3" << endl << nx << " " << ny << endl << 255 << endl;
    vector<Hitable *> list;
    list.push_back(new Sphere(Vector3d(0, 0, -1), 0.5, new Lambertian(Vector3d(0.8, 0.3, 0.3))));
    list.push_back(new Sphere(Vector3d(0, -100.5, -1), 100, new Lambertian(Vector3d(0.8, 0.8, 0))));
    list.push_back(new Sphere(Vector3d(1, 0, -1), 0.5, new Metal(Vector3d(0.8, 0.6, 0.2), 0.3)));
    list.push_back(new Sphere(Vector3d(-1, 0, -1), 0.5, new Dielectric(1.5)));
    Camera cam;
    auto world = new Hitable_list(list, list.size());
    for (int j = ny - 1; j >= 0; j--) {
        for (int i = 0; i < nx; i++) {
            Vector3d c(0, 0, 0);
            for (int s = 0; s < ns; s++) {
                double u = double(i + double(rand() % RAND_MAX) / double(RAND_MAX)) / double(nx);
                double v = double(j + double(rand() % RAND_MAX) / double(RAND_MAX)) / double(ny);
                Ray r = cam.get_ray(u, v);
//                Vector3d p = r.point_at_parameter(2.0);
                c += color(r, world, 0);
            }
            c /= double(ns);
            int ir = int(259.99 * sqrt(c(0)));
            int ig = int(259.99 * sqrt(c(1)));
            int ib = int(259.99 * sqrt(c(2)));
            outputFile << ir << " " << ig << " " << ib << endl;
        }
    }
    outputFile.close();
}