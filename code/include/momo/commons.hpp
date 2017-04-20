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
namespace commons {

template <typename Type>
static Eigen::Matrix<Type, 3, 3> GetCrossProdMatrix(Eigen::Matrix<Type, 3, 1> T) {
    Eigen::Matrix<Type, 3, 3> Out;
    Out << Type(0.), -T(2, 0), T(1, 0), T(2, 0), Type(0.), -T(0, 0), -T(1, 0), T(0, 0), Type(0.);

    return Out;
}

template <typename T>
T Squ(const T& a) {
    return a * a;
}

template <typename T>
T NormSquared(const std::pair<T, T>& a, const std::pair<T, T>& b) {
    return squ(a.first - b.first) + squ(a.second - b.second);
}

template <typename T>
T Norm(const std::pair<T, T>& a, const std::pair<T, T>& b) {
    return std::sqrt(norm_squared(a, b));
}

// scoped enum for telling the compiler which cost function to use
// for symmetric epipolar distance you need to assign seperate residual blocks with
// Forward("Normal") and Backward Epipolar Distance(all image domains)
// RayToEpipolarPlane is sin(a) where "a" is the deviation angle of the current vieweing ray to
// the epipolar plane(world domaine)
enum class CostFunctionType {
    SampsonDistance,
    EpipolarDistanceForward,
    EpipolarDistanceBackward,
    RayToEpipolarPlaneLinear,
    RayToEpipolarPlaneQuadratic
};
}
}
