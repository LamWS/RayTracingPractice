//
// Created by LamWS on 2021/5/14.
//

#ifndef RT_HITABLE_H
#define RT_HITABLE_H

#include <utility>

#include "Ray.h"
#include "Aabb.h"
#include "Eigen/Dense"
#include "algorithm"
#include <vector>

class Material;

struct hit_record {
    double t;
    Eigen::Vector3d p, normal;
    Material *material;
};

class Hitable {
public:
    virtual bool hit(const Ray &r, double t_min, double t_max, hit_record &rec) const = 0;

    virtual bool bounding_box(double t0, double t1, Aabb &box) const = 0;
};

class Sphere : public Hitable {
public:
    Sphere(Eigen::Vector3d cen, double r, Material *m) : center(std::move(cen)), material(m), radius(r) {}

    bool hit(const Ray &r, double t_min, double t_max, hit_record &rec) const override;

    bool bounding_box(double t0, double t1, Aabb &box) const override;

private:
    Eigen::Vector3d center;
    Material *material;
    double radius;
};

class Hitable_list : public Hitable {
public:
    Hitable_list(std::vector<Hitable *> l, int n);

    bool hit(const Ray &r, double t_min, double t_max, hit_record &rec) const override;

    bool bounding_box(double t0, double t1, Aabb &box) const override;

private:
    std::vector<Hitable *> list;
    int list_size;
};

class MovingSphere : public Hitable {
public:
//    MovingSphere() = default;

    MovingSphere(Eigen::Vector3d cen0, Eigen::Vector3d cen1, double t0, double t1, double r, Material *m) :
            center0(std::move(cen0)),
            center1(std::move(cen1)),
            time0(t0),
            time1(t1),
            radius(r),
            material(m) {}

    bool hit(const Ray &r, double t_min, double t_max, hit_record &rec) const override;

    bool bounding_box(double t0, double t1, Aabb &box) const override;

    Eigen::Vector3d center(double t) const;

private:
    Eigen::Vector3d center0, center1;
    double time0, time1;
    double radius;
    Material *material;
};

class BVHNode : public Hitable {
public:
    BVHNode() = default;

    BVHNode(std::vector<Hitable *> list, double t0, double t1);

    bool hit(const Ray &r, double t_min, double t_max, hit_record &rec) const override;

    bool bounding_box(double t0, double t1, Aabb &b) const override;

private:
    Hitable *left;
    Hitable *right;
    Aabb box;
};


#endif //RT_HITABLE_H
