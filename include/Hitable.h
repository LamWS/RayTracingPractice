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
    double u, v;
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

class XYRect : public Hitable {
public:
//    XYRect();
    XYRect(double x0, double x1, double y0, double y1, double k, Material *m) : x0(x0), x1(x1), y0(y0), y1(y1), k(k),
                                                                                material(m) {}

    bool hit(const Ray &r, double t_min, double t_max, hit_record &rec) const override;

    bool bounding_box(double t0, double t1, Aabb &box) const override;

private:
    double x0, x1, y0, y1, k;
    Material *material;
};

class XZRect : public Hitable {
public:
//    XYRect();
    XZRect(double x0, double x1, double z0, double z1, double k, Material *m) : x0(x0), x1(x1), z0(z0), z1(z1), k(k),
                                                                                material(m) {}

    bool hit(const Ray &r, double t_min, double t_max, hit_record &rec) const override;

    bool bounding_box(double t0, double t1, Aabb &box) const override;

private:
    double x0, x1, z0, z1, k;
    Material *material;
};

class YZRect : public Hitable {
public:
//    XYRect();
    YZRect(double x0, double x1, double z0, double z1, double k, Material *m) : y0(x0), y1(x1), z0(z0), z1(z1), k(k),
                                                                                material(m) {}

    bool hit(const Ray &r, double t_min, double t_max, hit_record &rec) const override;

    bool bounding_box(double t0, double t1, Aabb &box) const override;

private:
    double y0, y1, z0, z1, k;
    Material *material;
};

class FlipNormal : public Hitable {
public:
    explicit FlipNormal(Hitable *h) : hitable(h) {}

    bool hit(const Ray &r, double t_min, double t_max, hit_record &rec) const override;

    bool bounding_box(double t0, double t1, Aabb &box) const override;

    Hitable *hitable;
};

class Box : public Hitable {
public:
    Box(Eigen::Vector3d p0, Eigen::Vector3d p1, Material *m);

    bool hit(const Ray &r, double t_min, double t_max, hit_record &rec) const override;

    bool bounding_box(double t0, double t1, Aabb &box) const override;

private:
    Eigen::Vector3d p_min, p_max;
    Hitable *list;
};

#endif //RT_HITABLE_H
