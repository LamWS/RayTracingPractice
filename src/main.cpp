#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include "Ray.h"
#include "Hitable.h"
#include <cfloat>
#include "Camera.h"
#include "random"
#include "Material.h"
#include "omp.h"

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

double ran() {
    return ((rand() % RAND_MAX) / (double) RAND_MAX);
}

Hitable *random_scene() {
    int n = 50000;
    vector<Hitable *> list;
    list.push_back(new Sphere(Vector3d(0, -1000, 0), 1000, new Lambertian(Vector3d(0.5, 0.5, 0.5))));
    for (int a = -11; a < 11; a++) {
        for (int b = -11; b < 11; b++) {
            Vector3d center(a + 0.9 * ran(), 0.2,
                            b + 0.9 * ran());
            double choose_mat = ran();
            if ((center - Vector3d(4, 0.2, 0)).norm() > 0.9) {
                if (choose_mat < 0.8) {
//                    list.push_back(new Sphere(center, 0.2,
//                                              new Lambertian(
//                                                      Vector3d(0.5 * (1 + ran()), 0.5 * (1 + ran()), 0.5 * (1 + ran())))));
                    double t = 0.5 * ran();
                    Vector3d add(0, t, 0);
                    auto center1 = center + add;
                    cout << t << endl;
                    cout << "cen0" << endl;
                    cout << center << endl;
                    cout << "cen1" << endl;
                    cout << center1 << endl;
                    cout << "norm" << endl;
                    cout << (center1 - center).norm() << endl;
                    list.push_back(new MovingSphere(center, center1, 0, 1, 0.2,
                                                    new Lambertian(
                                                            Vector3d(ran() * ran(), ran() * ran(), ran() * ran()))));
                } else if (choose_mat < 0.95) {
                    list.push_back(new Sphere(center, 0.2,
                                              new Metal(
                                                      Vector3d(0.5 * (1 + ran()), 0.5 * (1 + ran()), 0.5 * (1 + ran())),
                                                      0.5 * ran())));
                } else {
                    list.push_back(new Sphere(center, 0.2, new Dielectric(1.5)));
                }
            }
        }
    }
    list.push_back(new Sphere(Vector3d(0, 1, 0), 1.0, new Dielectric(1.5)));
    list.push_back(new Sphere(Vector3d(-4, 1, 0), 1.0, new Lambertian(Vector3d(0.4, 0.2, 0.1))));
    list.push_back(new Sphere(Vector3d(4, 1, 0), 1.0, new Metal(Vector3d(0.7, 0.6, 0.5), 0)));
    return new Hitable_list(list, list.size());
}

int main() {
    ofstream outputFile;
    outputFile.open("output.ppm", ios::out);

    int nx = 200, ny = 100, ns = 100;
    outputFile << "P3" << endl << nx << " " << ny << endl << 255 << endl;
    vector<Hitable *> list;
    list.push_back(new Sphere(Vector3d(0, 0, -1), 0.5, new Lambertian(Vector3d(0.1, 0.2, 0.5))));
    list.push_back(new Sphere(Vector3d(0, -100.5, -1), 100, new Lambertian(Vector3d(0.8, 0.8, 0))));
    list.push_back(new Sphere(Vector3d(1, 0, -1), 0.5, new Metal(Vector3d(0.8, 0.6, 0.2), 0.3)));
    list.push_back(new Sphere(Vector3d(-1, 0, -1), 0.5, new Dielectric(1.5)));
    list.push_back(new Sphere(Vector3d(-1, 0, -1), -0.45, new Dielectric(1.5)));
//    double R = cos(M_PI / 4);
//    list.push_back(new Sphere(Vector3d(-R, 0, -1), R, new Lambertian(Vector3d(0, 0, 1))));
//    list.push_back(new Sphere(Vector3d(R, 0, -1), R, new Lambertian(Vector3d(1, 0, 0))));
    Vector3d lookFrom(13, 2, 3), lookAt(0, 0, 0);
    double dist_to_focus = (lookFrom - lookAt).norm();
    double aperture = 0;
    Camera cam(lookFrom, lookAt, Vector3d(0, 1, 0), 20, double(nx) / double(ny), aperture, dist_to_focus, 0, 1);
//    auto world = random_scene();
    auto world = new Hitable_list(list, list.size());
    vector<Vector3d> result;
    result.resize(ny * nx);
//#pragma omp parallel for default(none) shared(ny, nx, ns, world, cam, result)
    for (int j = ny - 1; j >= 0; j--) {
        cout << j << endl;
        for (int i = 0; i < nx; i++) {
//            cout << ny - j << ":" << i << endl;
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
//            result[(ny - j - 1) * nx + i] = c;
//            cout << (ny - j - 1) * nx + i << endl;
            outputFile << ir << " " << ig << " " << ib << endl;
        }
    }
//    for (auto c:result) {
//        int ir = int(259.99 * sqrt(c(0)));
//        int ig = int(259.99 * sqrt(c(1)));
//        int ib = int(259.99 * sqrt(c(2)));
//        outputFile << ir << " " << ig << " " << ib << endl;
//    }
    outputFile.close();
}