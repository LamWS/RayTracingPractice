//
// Created by LamWS on 2021/5/14.
//

#ifndef RT_CAMERA_H
#define RT_CAMERA_H

#include <Eigen/Dense>
#include "Ray.h"

class Camera {
public:
    Camera(Eigen::Vector3d lookFrom, const Eigen::Vector3d &lookAt, const Eigen::Vector3d &vup, double vfov,
           double aspect,
           double aperture, double focus_dist, double t0, double t1);

    Ray get_ray(double u, double v);

    Ray get_ray(Eigen::Vector3d look_at);

    Eigen::Vector3d screen_pos(const Eigen::Vector3d& look_at);

private:
    Eigen::Vector3d origin, lower_left_corner, horizontal, vertical, u, v, w;
    double len_radius;
    double time0, time1;
};


#endif //RT_CAMERA_H
