//
// Created by LamWS on 2021/5/13.
//

#ifndef RT_RAY_H
#define RT_RAY_H

#include <Eigen/Dense>

class Ray {
public:
    Ray(Eigen::Vector3d a, Eigen::Vector3d b);

    Eigen::Vector3d origin() const;

    Eigen::Vector3d direction() const;

    Eigen::Vector3d point_at_parameter(double t) const;

private:
    Eigen::Vector3d A, B;
};


#endif //RT_RAY_H
