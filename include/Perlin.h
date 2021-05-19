//
// Created by LamWS on 2021/5/19.
//

#ifndef RT_PERLIN_H
#define RT_PERLIN_H

#include <Eigen/Dense>

class Perlin {
public:
    double noise(const Eigen::Vector3d &p) const;

private:
    static Eigen::Vector3d *randVec;
    static int *perm_x;
    static int *perm_y;
    static int *perm_z;
};


#endif //RT_PERLIN_H
