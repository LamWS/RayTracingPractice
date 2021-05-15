//
// Created by LamWS on 2021/5/14.
//

#ifndef RT_CAMERA_H
#define RT_CAMERA_H

#include <Eigen/Dense>
#include "Ray.h"

class Camera {
public:
    Camera();

    Ray get_ray(double u, double v);

private:
    Eigen::Vector3d origin, lower_left_corner, horizontal, vertical;
};


#endif //RT_CAMERA_H
