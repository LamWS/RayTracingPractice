//
// Created by LamWS on 2021/5/14.
//

#include "Camera.h"

Camera::Camera(Eigen::Vector3d lookFrom, Eigen::Vector3d lookAt, Eigen::Vector3d vup, double vfov, double aspect) {
    Eigen::Vector3d u, v, w;
    double theta = vfov * M_PI / 180;
    double half_height = tan(theta / 2);
    double half_width = aspect * half_height;
    origin = std::move(lookFrom);
    w = (lookFrom - lookAt).normalized();
    u = vup.cross(w).normalized();
    v = w.cross(u);
//    lower_left_corner = Eigen::Vector3d(-half_width, -half_height, -1.0);
    lower_left_corner = origin - half_width * u - half_height * v - w;
    horizontal = 2 * half_width * u;
    vertical = 2 * half_height * v;
}

Ray Camera::get_ray(double u, double v) {
    return Ray(origin, lower_left_corner + u * horizontal + v * vertical - origin);
}
