//
// Created by lamws on 2021/5/19.
//

#ifndef RT_AABB_H
#define RT_AABB_H

#include <Eigen/Dense>
#include <utility>
#include "Ray.h"

class Aabb {
public:
    Aabb() = default;

    Aabb(Eigen::Vector3d maxx, Eigen::Vector3d minn) : _max(std::move(maxx)), _min(std::move(minn)) {}

    Eigen::Vector3d get_min() const;

    Eigen::Vector3d get_max() const;

    bool hit(const Ray &ray, double tmin, double tmax) const;

private:
    Eigen::Vector3d _min, _max;
};


#endif //RT_AABB_H
