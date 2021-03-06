//
// Created by LamWS on 2021/5/15.
//

#ifndef RT_MATERIAL_H
#define RT_MATERIAL_H

#include <utility>

#include "Ray.h"
#include "Hitable.h"
#include "Eigen/Dense"
#include "Texture.h"

class Material {
public:
    virtual bool scatter(const Ray &r_in, const hit_record &rec, Eigen::Vector3d &attenuation, Ray &scatter) const = 0;
    virtual Eigen::Vector3d emitted(double u, double v, const Eigen::Vector3d &p) const;
};

class Lambertian : public Material {
public:
    explicit Lambertian(Texture *a);

    bool scatter(const Ray &r_in, const hit_record &rec, Eigen::Vector3d &attenuation, Ray &scatter) const override;

private:
    Texture *albedo;
};

class Metal : public Material {
public:
    Metal(Eigen::Vector3d a, double f) : albedo(std::move(a)), fuzz(f) {}

    bool scatter(const Ray &r_in, const hit_record &rec, Eigen::Vector3d &attenuation, Ray &scatter) const override;

private:
    Eigen::Vector3d albedo;
    double fuzz;
};

class Dielectric : public Material {
public:
    explicit Dielectric(double ri) : ref_idx(ri) {}

    bool scatter(const Ray &r_in, const hit_record &rec, Eigen::Vector3d &attenuation, Ray &scatter) const override;

private:
    double ref_idx;
};

class DiffuseLight : public Material {
public:
    DiffuseLight(Texture *a) : emit(a) {}

    bool scatter(const Ray &r_in, const hit_record &rec, Eigen::Vector3d &attenuation, Ray &scatter) const override;

    Eigen::Vector3d emitted(double u, double v, const Eigen::Vector3d &p) const override;

private:
    Texture *emit;
};

#endif //RT_MATERIAL_H
