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

double schlick(double cosine, double ref_idx) {
    double r0 = (1 - ref_idx) / (1 + ref_idx);
    r0 = r0 * r0;
    return r0 + (1 - r0) * pow((1 - cosine), 5);
}

bool refract(const Vector3d &v, const Vector3d &n, double ni_over_nt, Vector3d &refracted) {
    Vector3d uv = v.normalized();
    double dt = uv.dot(n);
    double d = 1.0 - ni_over_nt * ni_over_nt * (1 - dt * dt);
    if (d > 0) {
        refracted = ni_over_nt * (uv - n * dt) - n * sqrt(d);
        return true;
    }
    return false;
}

bool Dielectric::scatter(const Ray &r_in, const hit_record &rec, Vector3d &attenuation, Ray &scatter) const {
    Vector3d outward_normal;
    Vector3d reflected = reflect(r_in.direction(), rec.normal);
    double ni_over_nt;
    attenuation = Vector3d(1.0, 1.0, 1.0);
    Vector3d refracted;

    // handle reflection using Schlick's approximation

    double reflect_prob;
    double cosine;

    if (r_in.direction().dot(rec.normal) > 0) {
        // surface downside in
        outward_normal = -rec.normal;
        ni_over_nt = ref_idx;
        cosine = ref_idx * r_in.direction().dot(rec.normal) / r_in.direction().squaredNorm();
    } else {
        outward_normal = rec.normal;
        ni_over_nt = 1 / ref_idx;
        cosine = -r_in.direction().dot(rec.normal) / r_in.direction().squaredNorm();
    }
    if (refract(r_in.direction(), outward_normal, ni_over_nt, refracted)) {
        reflect_prob = schlick(cosine, ref_idx);
    } else {
        reflect_prob = 1.0;
    }
    if ((rand() % RAND_MAX) / (double) RAND_MAX < reflect_prob) {
        scatter = Ray(rec.p, reflected);
    } else {
        scatter = Ray(rec.p, refracted);
    }
    return true;
}
