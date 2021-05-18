//
// Created by LamWS on 2021/5/13.
//

#include "Ray.h"

#include <utility>

Ray::Ray(Eigen::Vector3d a, Eigen::Vector3d b, double t) : A(std::move(a)), B(std::move(b)), _time(t) {}

Eigen::Vector3d Ray::origin() const {
    return A;
}

Eigen::Vector3d Ray::direction() const {
    return B;
}

Eigen::Vector3d Ray::point_at_parameter(double t) const {
    return A + t * B;
}

double Ray::time() const {
    return _time;
}

Ray::Ray() = default;
