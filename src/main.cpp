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

Hitable *cornell_box() {
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
    list.push_back(new Box(Vector3d(130, 0, 65), Vector3d(295, 165, 230), white));
    list.push_back(new Box(Vector3d(265, 0, 295), Vector3d(430, 330, 460), white));
    return new Hitable_list(list, list.size());
}

int main() {
    ofstream outputFile;
    outputFile.open("output.ppm", ios::out);

    int nx = 256, ny = 256, ns = 200;
    outputFile << "P3" << endl << nx << " " << ny << endl << 255 << endl;
//    vector<Hitable *> list;
//    list.push_back(new Sphere(Vector3d(0, 0, -1), 1, new Lambertian(new ConstantTexture(Vector3d(0.1, 0.2, 0.5)))));
//    list.push_back(
//            new Sphere(Vector3d(0, -100.5, -1), 100, new Lambertian(new ConstantTexture(Vector3d(0.8, 0.8, 0)))));
//    list.push_back(new Sphere(Vector3d(1, 0, -1), 1, new Metal(Vector3d(0.8, 0.6, 0.2), 0.3)));
//    list.push_back(new Sphere(Vector3d(-1, 0, -1), 1, new Dielectric(1.5)));
//    list.push_back(new Sphere(Vector3d(0, 7, 0), 2, new DiffuseLight(new ConstantTexture(Vector3d(4, 4, 4)))));
//    list.push_back(new XYRect(3, 5, 1, 3, -2, new DiffuseLight(new ConstantTexture(Vector3d(4, 4, 4)))));
//    list.push_back(new Sphere(Vector3d(-1, 0, -1), -0.45, new Dielectric(1.5)));
//    double R = cos(M_PI / 4);
//    list.push_back(new Sphere(Vector3d(-R, 0, -1), R, new Lambertian(Vector3d(0, 0, 1))));
//    list.push_back(new Sphere(Vector3d(R, 0, -1), R, new Lambertian(Vector3d(1, 0, 0))));
//    list.push_back(new Sphere(Vector3d(0, -1000, 0), 1000, new Lambertian(new NoiseTexture(1))));
//    list.push_back(new Sphere(Vector3d(0, 2, 0), 2, new Lambertian(new NoiseTexture(1))));
    Vector3d lookFrom(278, 278, -800), lookAt(278, 278, 0);
    double dist_to_focus = 10;
    double aperture = 0;
    Camera cam(lookFrom, lookAt, Vector3d(0, 1, 0), 40, double(nx) / double(ny), aperture, dist_to_focus, 0, 1);
    auto world = cornell_box();
//    auto world = new Hitable_list(list, list.size());
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