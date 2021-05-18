//
// Created by LamWS on 2021/5/13.
//

#ifndef RT_RAY_H
#define RT_RAY_H

#include <Eigen/Dense>

class Ray {
public:
    Ray();
    Ray(Eigen::Vector3d a, Eigen::Vector3d b,double t);

    Eigen::Vector3d origin() const;

    Eigen::Vector3d direction() const;
    double time() const;
    Eigen::Vector3d point_at_parameter(double t) const;

private:
    Eigen::Vector3d A, B;
    double _time;
};


#endif //RT_RAY_H
