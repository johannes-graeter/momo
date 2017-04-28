/*
 * Copyright 2017. All rights reserved.
 * Institute of Measurement and Control Systems
 * Karlsruhe Institute of Technology, Germany
 *
 * authors:
 *  Johannes Graeter (johannes.graeter@kit.edu)
 *  and others
 */
#pragma once
#include <Eigen/Eigen>

namespace momo {
namespace ceres_eigen_types {

template<typename T>
using Affine3 = Eigen::Transform<T, 3, Eigen::Affine, Eigen::DontAlign>;

template<typename T>
using Matrix3 = Eigen::Matrix<T, 3, 3, Eigen::DontAlign>;

template<typename T>
using Vector3 = Eigen::Matrix<T, 3, 1, Eigen::DontAlign>;
}
}
