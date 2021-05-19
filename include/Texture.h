//
// Created by LamWS on 2021/5/19.
//

#ifndef RT_TEXTURE_H
#define RT_TEXTURE_H

#include <Eigen/Dense>
#include <utility>

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

#endif //RT_TEXTURE_H
