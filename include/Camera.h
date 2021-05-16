//
// Created by LamWS on 2021/5/14.
//

#ifndef RT_CAMERA_H
#define RT_CAMERA_H

#include <Eigen/Dense>
#include "Ray.h"

class Camera {
public:
    Camera(Eigen::Vector3d lookFrom, Eigen::Vector3d lookAt, Eigen::Vector3d vup, double vfov, double aspect,
           double aperture, double focus_dist);

    Ray get_ray(double u, double v);

private:
    Eigen::Vector3d origin, lower_left_corner, horizontal, vertical, u, v, w;
    double len_radius;

};


#endif //RT_CAMERA_H
