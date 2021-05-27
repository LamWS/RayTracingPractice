//
// Created by LamWS on 2021/5/26.
//

#ifndef RT_EDGE_H
#define RT_EDGE_H

#include <Eigen/Dense>
#include <utility>

class Edge {
public:
    Edge(Eigen::Vector3d v0, Eigen::Vector3d v1) : v0(std::move(v0)), v1(std::move(v1)) {}

    Eigen::Vector3d v0, v1;
};


#endif //RT_EDGE_H
