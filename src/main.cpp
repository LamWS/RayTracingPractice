#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include "Ray.h"
#include "Hitable.h"
#include <cfloat>
#include "Camera.h"
#include "random"

using namespace std;
using namespace Eigen;

Vector3d random_in_unit_sphere() {
    Vector3d p;
    do {
        double a = double(rand() % RAND_MAX) / double(RAND_MAX), b =
                double(rand() % RAND_MAX) / double(RAND_MAX), c = double(rand() % RAND_MAX) / double(RAND_MAX);
        p = 2.0 * Vector3d(a, b, c) - Vector3d(1, 1, 1);
    } while (p.x() * p.x() + p.y() * p.y() + p.z() * p.z() >= 1);
    return p;
}

Vector3d color(const Ray &r, Hitable *world) {
    hit_record rec;
    if (world->hit(r, 0.001, FLT_MAX, rec)) {
        Vector3d target = rec.p + rec.normal + random_in_unit_sphere();
        return 0.5 * color(Ray(rec.p, target - rec.p), world);
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
    list.push_back(new Sphere(Vector3d(0, 0, -1), 0.5));
    list.push_back(new Sphere(Vector3d(0, -100.5, -1), 100));
    Camera cam;
    auto world = new Hitable_list(list, 2);
    for (int j = ny - 1; j >= 0; j--) {
        for (int i = 0; i < nx; i++) {
            Vector3d c(0, 0, 0);
            for (int s = 0; s < ns; s++) {
                double u = double(i + double(rand() % RAND_MAX) / double(RAND_MAX)) / double(nx);
                double v = double(j + double(rand() % RAND_MAX) / double(RAND_MAX)) / double(ny);
                Ray r = cam.get_ray(u, v);
//                Vector3d p = r.point_at_parameter(2.0);
                c += color(r, world);
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