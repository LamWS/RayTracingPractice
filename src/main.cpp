#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include "Ray.h"
#include "Hitable.h"
#include <cfloat>
#include <Edge.h>
#include "Camera.h"
#include "random"
#include "Material.h"
#include "omp.h"
#include "Texture.h"

#define STB_IMAGE_IMPLEMENTATION

#include "stb_image.h"

using namespace std;
using namespace Eigen;


Vector3d color(const Ray &r, Hitable *world, int depth) {
    hit_record rec;
    if (world->hit(r, 0.001, FLT_MAX, rec)) {
        Ray scatter;
        Vector3d attenuation;
        auto emitted = rec.material->emitted(rec.u, rec.v, rec.p);
        if (depth < 50 && rec.material->scatter(r, rec, attenuation, scatter)) {
            auto c = color(scatter, world, depth + 1);
            return emitted + Vector3d(c.x() * attenuation.x(), c.y() * attenuation.y(), c.z() * attenuation.z());
        }
        return emitted;
    }
    return Vector3d(0, 0, 0);
}

double ran() {
    return ((rand() % RAND_MAX) / (double) RAND_MAX);
}

Hitable *random_scene() {
    int n = 50000;
    vector<Hitable *> list;
    list.push_back(
            new Sphere(Vector3d(0, -1000, 0), 1000, new Lambertian(
                    new CheckerTexture(new ConstantTexture(Vector3d(0.2, 0.3, 0.1)),
                                       new ConstantTexture(Vector3d(0.9, 0.9, 0.9))
                    ))
            )
    );
    for (int a = -11; a < 11; a++) {
        for (int b = -11; b < 11; b++) {
            Vector3d center(a + 0.9 * ran(), 0.2,
                            b + 0.9 * ran());
            double choose_mat = ran();
            if ((center - Vector3d(4, 0.2, 0)).norm() > 0.9) {
                if (choose_mat < 0.8) {
                    double t = 0.5 * ran();
                    Vector3d add(0, t, 0);
                    auto center1 = center + add;
                    list.push_back(new MovingSphere(center, center1, 0, 1, 0.2,
                                                    new Lambertian(new ConstantTexture(
                                                            Vector3d(ran() * ran(),
                                                                     ran() * ran(),
                                                                     ran() * ran())
                                                    ))));
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
    list.push_back(new Sphere(Vector3d(-4, 1, 0), 1.0, new Lambertian(new ConstantTexture(Vector3d(0.4, 0.2, 0.1)))));
    list.push_back(new Sphere(Vector3d(4, 1, 0), 1.0, new Metal(Vector3d(0.7, 0.6, 0.5), 0)));
    return new Hitable_list(list, list.size());
}

Hitable *cornell_box(double eps = 0) {
    vector<Hitable *> list;
    auto red = new Lambertian(new ConstantTexture(Vector3d(0.65, 0.05, 0.05)));
    auto white = new Lambertian(new ConstantTexture(Vector3d(0.73, 0.73, 0.73)));
    auto green = new Lambertian(new ConstantTexture(Vector3d(0.12, 0.45, 0.15)));
    auto light = new DiffuseLight(new ConstantTexture(Vector3d(15, 15, 15)));
    list.push_back(new FlipNormal(new YZRect(0, 555, 0, 555, 555, green)));
    list.push_back(new YZRect(0, 555, 0, 555, 0, red));
    list.push_back(new XZRect(213, 343, 227, 332, 554, light));
    list.push_back(new FlipNormal(new XZRect(0, 555, 0, 555, 555, white)));
    list.push_back(new XZRect(0, 555, 0, 555, 0, white));
    list.push_back(new FlipNormal(new XYRect(0, 555, 0, 555, 555, white)));
//    list.push_back(new Box(Vector3d(130, 0, 65), Vector3d(295, 165, 230), white));
//    list.push_back(new Box(Vector3d(265 + eps, 0, 150), Vector3d(430 + eps, 200, 380), white));
//    list.push_back(new Box(Vector3d(130, 0, 65), Vector3d(295, 165, 230), white));
//    list.push_back(new Box(Vector3d(265, 0, 295), Vector3d(430, 330, 460), white));
    list.push_back(new Sphere(Vector3d(200, 200, 200), 100.0, new Metal(Vector3d(0.7, 0.6, 0.5), 0)));
    return new Hitable_list(list, list.size());
}

int main() {
    ofstream outputFile;
    outputFile.open("output.ppm", ios::out);
    int nx = 256, ny = 256, ns = 6400;
    outputFile << "P3" << endl << nx << " " << ny << endl << 255 << endl;
    Vector3d lookFrom(278, 278, -800), lookAt(278, 278, 0);
    Vector3d x_eps(0.01, 0, 0);
    double dist_to_focus = 10;
    double aperture = 0;
    Camera cam(lookFrom, lookAt, Vector3d(0, 1, 0), 40, double(nx) / double(ny), aperture, dist_to_focus, 0, 1);
    Camera cam1(lookFrom + x_eps, lookAt + x_eps, Vector3d(0, 1, 0), 40, double(nx) / double(ny), aperture,
                dist_to_focus, 0, 1);
    auto world = cornell_box();
    vector<Vector3d> result;
    result.resize(ny * nx);
    for (int j = ny - 1; j >= 0; j--) {
        cout << j << endl;
        for (int i = 0; i < nx; i++) {
            Vector3d c(0, 0, 0);
            for (int s = 0; s < ns; s++) {
                double u = double(i + double(rand() % RAND_MAX) / double(RAND_MAX)) / double(nx);
                double v = double(j + double(rand() % RAND_MAX) / double(RAND_MAX)) / double(ny);
                Ray r = cam.get_ray(u, v);
                c += color(r, world, 0);
            }
//            for (int s = 0; s < ns; s++) {
//                double u = double(i + double(rand() % RAND_MAX) / double(RAND_MAX)) / double(nx);
//                double v = double(j + double(rand() % RAND_MAX) / double(RAND_MAX)) / double(ny);
//                Ray r = cam1.get_ray(u, v);
//                c_out += color(r, world, 0);
//            }
            c /= double(ns);
//            c_out /= double(ns);
//            auto c = (c_out - c_in) / x_eps.x();
            int ir = int(259.99 * sqrt(clamp(c(0), 0.0, 1.0)));
            int ig = int(259.99 * sqrt(clamp(c(1), 0.0, 1.0)));
            int ib = int(259.99 * sqrt(clamp(c(2), 0.0, 1.0)));
            outputFile << ir << " " << ig << " " << ib << endl;
        }
    }
    outputFile.close();
//    vector<Edge> edges;
//    edges.emplace_back(Vector3d(130, 0, 230), Vector3d(295, 0, 230));
//    edges.emplace_back(Vector3d(130, 0, 230), Vector3d(130, 165, 230));
//    edges.emplace_back(Vector3d(130, 165, 230), Vector3d(295, 165, 230));
//    edges.emplace_back(Vector3d(295, 0, 230), Vector3d(295, 165, 230));
//
//    edges.emplace_back(Vector3d(130, 0, 65), Vector3d(130, 0, 230));
//    edges.emplace_back(Vector3d(130, 165, 65), Vector3d(130, 165, 230));
//    edges.emplace_back(Vector3d(295, 0, 65), Vector3d(295, 0, 230));
//    edges.emplace_back(Vector3d(295, 165, 65), Vector3d(295, 165, 230));
//
//
//    edges.emplace_back(Vector3d(130, 0, 65), Vector3d(295, 0, 65));
//    edges.emplace_back(Vector3d(130, 0, 65), Vector3d(130, 165, 65));
//    edges.emplace_back(Vector3d(130, 165, 65), Vector3d(295, 165, 65));
//    edges.emplace_back(Vector3d(295, 0, 65), Vector3d(295, 165, 65));
//
//
//
//    edges.emplace_back(Vector3d(265, 0, 460), Vector3d(430, 0, 460));
//    edges.emplace_back(Vector3d(265, 0, 460), Vector3d(265, 330, 460));
//    edges.emplace_back(Vector3d(265, 330, 460), Vector3d(430, 330, 460));
//    edges.emplace_back(Vector3d(430, 0, 460), Vector3d(430, 330, 460));
//
//    edges.emplace_back(Vector3d(265, 0, 295), Vector3d(265, 0, 460));
//    edges.emplace_back(Vector3d(265, 330, 295), Vector3d(265, 330, 460));
//    edges.emplace_back(Vector3d(430, 0, 295), Vector3d(430, 0, 460));
//    edges.emplace_back(Vector3d(430, 330, 295), Vector3d(430, 330, 460));
//
//
//    edges.emplace_back(Vector3d(265, 0, 295), Vector3d(430, 0, 295));
//    edges.emplace_back(Vector3d(265, 0, 295), Vector3d(265, 330, 295));
//    edges.emplace_back(Vector3d(265, 330, 295), Vector3d(430, 330, 295));
//    edges.emplace_back(Vector3d(430, 0, 295), Vector3d(430, 330, 295));
//
//
//    edges.emplace_back(Vector3d(0, 0, 555), Vector3d(0, 555, 555));
//    edges.emplace_back(Vector3d(0, 0, 555), Vector3d(555, 0, 555));
//    edges.emplace_back(Vector3d(555, 555, 555), Vector3d(0, 555, 555));
//    edges.emplace_back(Vector3d(555, 555, 555), Vector3d(555, 0, 555));
//
//    edges.emplace_back(Vector3d(0, 0, 0), Vector3d(0, 0, 555));
//    edges.emplace_back(Vector3d(0, 555, 0), Vector3d(0, 555, 555));
//    edges.emplace_back(Vector3d(555, 0, 0), Vector3d(555, 0, 555));
//    edges.emplace_back(Vector3d(555, 555, 0), Vector3d(555, 555, 555));
//
//
//    edges.emplace_back(Vector3d(0, 0, 0), Vector3d(0, 555, 0));
//    edges.emplace_back(Vector3d(0, 0, 0), Vector3d(555, 0, 0));
//    edges.emplace_back(Vector3d(555, 555, 0), Vector3d(0, 555, 0));
//    edges.emplace_back(Vector3d(555, 555, 0), Vector3d(555, 0, 0));
//
////    Vector3d(265 + eps, 0, 295), Vector3d(430 + eps, 330, 460)
//    vector<Vector3d> screen_dx, screen_dy, screen_dz;
//    screen_dx.resize(nx * ny);
//    screen_dy.resize(nx * ny);
//    screen_dz.resize(nx * ny);
//    int sample_num = 20000;
//    for (int i = 0; i < sample_num; i++) {
//        cout << i << endl;
//        int edge_id = rand() % edges.size();
//        auto edge = edges[edge_id];
//        auto v0 = edge.v0;
//        auto v1 = edge.v1;
//        auto t = ran();
////        cout << t << endl;
//        auto p = v0 + t * (v1 - v0);
//        auto direct = (v1 - v0).normalized();
//        Vector3d n(-direct.y(), direct.x(), direct.z());
////        cout << n << endl;
//        // t, u, v
//        auto pos = cam.screen_pos(p);
////        cout << pos << endl;
//        int ti = (int) pos.x(), xi = (int) (pos.y() * nx), yi = (int) (pos.z() * ny);
//        if (ti < 0 || ti > 1 || xi < 0 || xi >= nx || yi < 0 || yi >= ny)continue;
////        cout << p << endl;
//
//        Vector3d color_in(0, 0, 0), color_out(0, 0, 0);
//        for (int s = 0; s < ns; s++) {
//            Ray r_in = cam.get_ray(p - 1e-3f * n);
//            Ray r_out = cam.get_ray(p + 1e-3f * n);
//            color_in += color(r_in, world, 0);
//            color_out += color(r_out, world, 0);
//        }
//        color_in /= double(ns);
//        color_out /= double(ns);
//        auto pdf = 1 / (v1 - v0).norm();
//        auto weight = 1 / (pdf * (double) sample_num);
//        Vector3d dx = -n.x() * (color_in - color_out) * weight;
//        Vector3d dy = -n.y() * (color_in - color_out) * weight;
//        Vector3d dz = -n.z() * (color_in - color_out) * weight;
////        if ((color_in - color_out).norm() > 0.011)
////            cout << "cap" << endl;
////        cout << "in" << endl;
////        cout << color_in << endl;
////        cout << "out" << endl;
////        cout << color_in << endl;
////        auto dy = -n.y() * (color_in - color_out) * weight;
////        auto dz = -n.z() * (color_in - color_out) * weight;
////        cout << yi * nx + xi << endl;
////        cout << yi << endl;
//        screen_dx[(ny - yi - 1) * nx + xi] += dx;
//        screen_dy[(ny - yi - 1) * nx + xi] += dy;
//        screen_dz[(ny - yi - 1) * nx + xi] += dz;
//    }
//
//    outputFile.open("dx.ppm", ios::out);
//
//    outputFile << "P3" << endl << nx << " " << ny << endl << 255 << endl;
//    for (auto c:screen_dx) {
//        int ir = int(259.99 * sqrt(clamp(c(0), 0.0, 1.0)));
//        int ig = int(259.99 * sqrt(clamp(c(1), 0.0, 1.0)));
//        int ib = int(259.99 * sqrt(clamp(c(2), 0.0, 1.0)));
//        outputFile << ir << " " << ig << " " << ib << endl;
//    }
//    outputFile.close();
//    outputFile.open("dy.ppm", ios::out);
//
//    outputFile << "P3" << endl << nx << " " << ny << endl << 255 << endl;
//    for (auto c:screen_dy) {
//        int ir = int(259.99 * sqrt(clamp(c(0), 0.0, 1.0)));
//        int ig = int(259.99 * sqrt(clamp(c(1), 0.0, 1.0)));
//        int ib = int(259.99 * sqrt(clamp(c(2), 0.0, 1.0)));
//        outputFile << ir << " " << ig << " " << ib << endl;
//    }
//    outputFile.close();
//    outputFile.open("dz.ppm", ios::out);
//
//    outputFile << "P3" << endl << nx << " " << ny << endl << 255 << endl;
//    for (auto c:screen_dz) {
//        int ir = int(259.99 * sqrt(clamp(c(0), 0.0, 1.0)));
//        int ig = int(259.99 * sqrt(clamp(c(1), 0.0, 1.0)));
//        int ib = int(259.99 * sqrt(clamp(c(2), 0.0, 1.0)));
//        outputFile << ir << " " << ig << " " << ib << endl;
//    }
//    outputFile.close();
    return 0;
}