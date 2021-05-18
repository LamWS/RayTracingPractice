//
// Created by LamWS on 2021/5/14.
//

#include "Camera.h"

using namespace Eigen;

Vector3d random_in_unit_disk() {
    Vector3d p;
    do {
        double a = double(rand() % RAND_MAX) / double(RAND_MAX),
                b = double(rand() % RAND_MAX) / double(RAND_MAX);
        p = 2.0 * Vector3d(a, b, 0) - Vector3d(1, 1, 0);
    } while (p.norm() >= 1);
    return p;
}

Camera::Camera(Eigen::Vector3d lookFrom, const Eigen::Vector3d &lookAt, const Eigen::Vector3d &vup, double vfov,
               double aspect,
               double aperture, double focus_dist, double t0, double t1) {
    len_radius = aperture / 2;
    time0 = t0;
    time1 = t1;
    double theta = vfov * M_PI / 180;
    double half_height = tan(theta / 2);
    double half_width = aspect * half_height;
    origin = std::move(lookFrom);
    w = (origin - lookAt).normalized();
    u = vup.cross(w).normalized();
    v = w.cross(u);
    lower_left_corner = origin - half_width * focus_dist * u - half_height * focus_dist * v - w * focus_dist;
    horizontal = 2 * half_width * focus_dist * u;
    vertical = 2 * half_height * focus_dist * v;
}

Ray Camera::get_ray(double s, double t) {
    Vector3d rd = len_radius * random_in_unit_disk();
    Vector3d offset = u * rd.x() + v * rd.y();
    double time = time0 + ((rand() % RAND_MAX) / (double) RAND_MAX) * (time1 - time0);
    return Ray(origin + offset, lower_left_corner + s * horizontal + t * vertical - origin - offset, time);
}
