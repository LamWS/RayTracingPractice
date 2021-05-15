//
// Created by LamWS on 2021/5/15.
//

#ifndef RT_MATERIAL_H
#define RT_MATERIAL_H

#include <utility>

#include "Ray.h"
#include "Hitable.h"
#include "Eigen/Dense"

class Material {
public:
    virtual bool scatter(const Ray &r_in, const hit_record &rec, Eigen::Vector3d &attenuation, Ray &scatter) const = 0;
};

class Lambertian : public Material {
public:
    explicit Lambertian(Eigen::Vector3d a);

    bool scatter(const Ray &r_in, const hit_record &rec, Eigen::Vector3d &attenuation, Ray &scatter) const override;

private:
    Eigen::Vector3d albedo;
};

class Metal : public Material {
public:
    Metal(Eigen::Vector3d a, double f) : albedo(std::move(a)), fuzz(f) {}

    bool scatter(const Ray &r_in, const hit_record &rec, Eigen::Vector3d &attenuation, Ray &scatter) const override;

private:
    Eigen::Vector3d albedo;
    double fuzz;
};

#endif //RT_MATERIAL_H
