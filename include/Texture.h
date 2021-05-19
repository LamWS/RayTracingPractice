//
// Created by LamWS on 2021/5/19.
//

#ifndef RT_TEXTURE_H
#define RT_TEXTURE_H

#include <Eigen/Dense>
#include <utility>
#include "Perlin.h"
#include <algorithm>
class Texture {
public:
    virtual Eigen::Vector3d value(double u, double v, const Eigen::Vector3d &p) const = 0;
};

class ConstantTexture : public Texture {
public:
    ConstantTexture() = default;

    explicit ConstantTexture(Eigen::Vector3d c) : color(std::move(c)) {}

    Eigen::Vector3d value(double u, double v, const Eigen::Vector3d &p) const override;

private:
    Eigen::Vector3d color;
};

class CheckerTexture : public Texture {
public:
    CheckerTexture() = default;

    CheckerTexture(Texture *t0, Texture *t1) : even(t0), odd(t1) {}

    Eigen::Vector3d value(double u, double v, const Eigen::Vector3d &p) const override;

private:
    Texture *even;
    Texture *odd;
};

class NoiseTexture : public Texture {
public:
    NoiseTexture() = default;

    explicit NoiseTexture(double sc) : scale(sc) {}

    Eigen::Vector3d value(double u, double v, const Eigen::Vector3d &p) const override;

private:
    Perlin noise;
    double scale;
};

class ImageTexture : public Texture {
public:
    ImageTexture(unsigned char *pixels, int A, int B) : data(pixels), nx(A), ny(B) {}
    Eigen::Vector3d value(double u, double v, const Eigen::Vector3d &p) const override;
private:
    unsigned char *data;
    int nx, ny;
};

#endif //RT_TEXTURE_H
