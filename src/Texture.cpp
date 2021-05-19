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

Eigen::Vector3d ImageTexture::value(double u, double v, const Eigen::Vector3d &p) const {
    int i = std::min((int) std::max(0.0, u * nx), nx - 1);
    int j = std::min((int) std::max(0.0, (1 - v) * ny - 0.001), ny - 1);
    double r = int(data[3 * i + 3 * nx * j]) / 255.0;
    double g = int(data[3 * i + 3 * nx * j + 1]) / 255.0;
    double b = int(data[3 * i + 3 * nx * j + 2]) / 255.0;
    return Eigen::Vector3d(r, g, b);
}
