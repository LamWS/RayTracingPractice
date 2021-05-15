//
// Created by LamWS on 2021/5/14.
//

#ifndef RT_HITABLE_H
#define RT_HITABLE_H

#include <utility>

#include "Ray.h"
#include "Eigen/Dense"
#include <vector>

//class Material;

struct hit_record {
    double t;
    Eigen::Vector3d p, normal;

};

class Hitable {
public:
    virtual bool hit(const Ray &r, double t_min, double t_max, hit_record &rec) const = 0;
};

class Sphere : public Hitable {
public:
    Sphere(Eigen::Vector3d cen, double r) : center(std::move(cen)), radius(r) {}

    bool hit(const Ray &r, double t_min, double t_max, hit_record &rec) const override;

private:
    Eigen::Vector3d center;
    double radius;
};

class Hitable_list : public Hitable {
public:
    Hitable_list(std::vector<Hitable *> l, int n);
    bool hit(const Ray &r, double t_min, double t_max, hit_record &rec) const override;

private:
    std::vector<Hitable *> list;
    int list_size;
};

#endif //RT_HITABLE_H
