//
// Created by LamWS on 2021/5/19.
//

#include "Texture.h"

Eigen::Vector3d ConstantTexture::value(double u, double v, const Eigen::Vector3d &p) const {
    return color;
}

Eigen::Vector3d CheckerTexture::value(double u, double v, const Eigen::Vector3d &p) const {
    double sines = sin(10 * p.x()) * sin(10 * p.y()) * sin(10 * p.z());
    if (sines < 0) {
        return odd->value(u, v, p);
    } else {
        return even->value(u, v, p);
    }
}

Eigen::Vector3d NoiseTexture::value(double u, double v, const Eigen::Vector3d &p) const {
    return Eigen::Vector3d(1, 1, 1) * noise.noise(scale * p);
}
