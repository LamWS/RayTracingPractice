//
// Created by LamWS on 2021/5/14.
//

#include "Camera.h"

Camera::Camera() {
    lower_left_corner = Eigen::Vector3d(-2.0, -1.0, -1.0);
    horizontal = Eigen::Vector3d(4.0, 0.0, 0.0);
    vertical = Eigen::Vector3d(0.0, 2.0, 0.0);
    origin = Eigen::Vector3d(0.0, 0.0, 0.0);
}

Ray Camera::get_ray(double u, double v) {
    return Ray(origin, lower_left_corner + u * horizontal + v * vertical - origin);
}
