//
// Created by lamws on 2021/5/19.
//

#include "Aabb.h"

Eigen::Vector3d Aabb::get_min() const {
    return _min;
}

Eigen::Vector3d Aabb::get_max() const {
    return _max;
}

bool Aabb::hit(const Ray &ray, double tmin, double tmax) const {
    for (int dim = 0; dim < 3; dim++) {
        double t0 = fmin((_min(dim) - ray.origin()(dim)) / ray.direction()(dim),
                         (_max(dim) - ray.origin()(dim)) / ray.direction()(dim));
        double t1 = fmax((_min(dim) - ray.origin()(dim)) / ray.direction()(dim),
                         (_max(dim) - ray.origin()(dim)) / ray.direction()(dim));
        tmin = fmax(t0, tmin);
        tmax = fmax(t1, tmax);
        if (tmax <= tmin)return false;
    }
    return true;
}
