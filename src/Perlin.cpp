//
// Created by LamWS on 2021/5/19.
//

#include "Perlin.h"
#include <iostream>

static Eigen::Vector3d *perlin_generate() {
    auto p = new Eigen::Vector3d[256];
    for (int i = 0; i < 256; i++) {
        p[i] = Eigen::Vector3d(-1 + 2 * (rand() % RAND_MAX) / (double) RAND_MAX,
                               -1 + 2 * (rand() % RAND_MAX) / (double) RAND_MAX,
                               -1 + 2 * (rand() % RAND_MAX) / (double) RAND_MAX);
//        p[i] = (rand() % RAND_MAX) / (double) RAND_MAX;
    }
    return p;
}

void permute(int *p, int n) {
    for (int i = n - 1; i > 0; i--) {
        int target = rand() % (i + 1);
        int tmp = p[i];
        p[i] = p[target];
        p[target] = tmp;
    }
}

static int *perlin_generate_perm() {
    auto *p = new int[256];
    for (int i = 0; i < 256; i++) {
        p[i] = i;
    }
    permute(p, 256);
    return p;
}

inline double triLinear_interp(double c[2][2][2], double u, double v, double w) {
    double accum = 0;
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            for (int k = 0; k < 2; k++) {
                accum += (i * u + (1 - i) * (1 - u)) *
                         (j * v + (1 - j) * (1 - v)) *
                         (k * w + (1 - k) * (1 - w)) * c[i][j][k];
            }
        }
    }
    return accum;
}

inline double perlin_interp(Eigen::Vector3d c[2][2][2], double u, double v, double w) {
    double uu = u * u * (3 - 2 * u);
    double vv = v * v * (3 - 2 * v);
    double ww = w * w * (3 - 2 * w);
    double accum = 0;
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            for (int k = 0; k < 2; k++) {
                Eigen::Vector3d weight_v(u - i, v - j, w - k);
                accum += (i * uu + (1 - i) * (1 - uu)) *
                         (j * vv + (1 - j) * (1 - vv)) *
                         (k * ww + (1 - k) * (1 - ww)) * c[i][j][k].dot(weight_v);
            }
        }
    }
    return accum;
}

double Perlin::noise(const Eigen::Vector3d &p) const {
    double u = p.x() - floor(p.x());
    double v = p.y() - floor(p.y());
    double w = p.z() - floor(p.z());
//    u = u * u * (3 - 2 * u);
//    v = v * v * (3 - 2 * v);
//    w = w * w * (3 - 2 * w);
    int i = floor(p.x());
    int j = floor(p.y());
    int k = floor(p.z());
    Eigen::Vector3d c[2][2][2];
    for (int di = 0; di < 2; di++) {
        for (int dj = 0; dj < 2; dj++) {
            for (int dk = 0; dk < 2; dk++) {
                c[di][dj][dk] = randVec[perm_x[(i + di) & 255] ^ perm_y[(j + dj) & 255] ^ perm_z[(k + dk) & 255]];
            }
        }
    }

    return perlin_interp(c, u, v, w);
}


Eigen::Vector3d *Perlin::randVec = perlin_generate();
int *Perlin::perm_x = perlin_generate_perm();
int *Perlin::perm_y = perlin_generate_perm();
int *Perlin::perm_z = perlin_generate_perm();