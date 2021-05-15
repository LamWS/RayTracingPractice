//
// Created by LamWS on 2021/5/15.
//

#include "Material.h"

#include <utility>

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


bool Lambertian::scatter(const Ray &r_in, const hit_record &rec, Eigen::Vector3d &attenuation, Ray &scatter) const {
    Eigen::Vector3d target = rec.p + rec.normal + random_in_unit_sphere();
    scatter = Ray(rec.p, target - rec.p);
    attenuation = albedo;
    return true;
}

Lambertian::Lambertian(Eigen::Vector3d a) : albedo(std::move(a)) {

}

Vector3d reflect(const Vector3d &v, const Vector3d &n) {
    return v - 2 * v.dot(n) * n;
}

bool Metal::scatter(const Ray &r_in, const hit_record &rec, Vector3d &attenuation, Ray &scatter) const {
    Vector3d reflected = reflect(r_in.direction(), rec.normal).normalized();
    scatter = Ray(rec.p, reflected + fuzz * random_in_unit_sphere());
    attenuation = albedo;
    return scatter.direction().dot(rec.normal) > 0;
}
