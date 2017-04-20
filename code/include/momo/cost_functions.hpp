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
#include <ceres/jet.h>

namespace momo {

namespace cost_functions {

template <typename T>
T squ(const T& a) {
    return a * a;
}

///@class interface class for cost_funcros must be templated -> we can emulate that with the
/// curiously recurring template pattern
/// http://stackoverflow.com/questions/8611318/any-way-to-have-a-template-function-in-an-abstract-base-class
// class CostFunctorInterface{

//};

namespace sampson_distance {

template <typename T>
T get_error(const Eigen::Matrix<T, 3, 3, Eigen::DontAlign>& fundamental_matrix,
            const Eigen::Matrix<T, 3, 1, Eigen::DontAlign>& x_cur,
            const Eigen::Matrix<T, 3, 1, Eigen::DontAlign>& x_old) {
    // calc sampson dist(Mulptiple view geom p.287)

    // denominator
    Eigen::Matrix<T, 3, 1, Eigen::DontAlign> Fx = fundamental_matrix * x_old;
    // Fx is epipolar line, if its norm is too small we will have an unstable error since it will be
    // normed afterwards -> return a penalty
    ///@todo
    Eigen::Matrix<T, 3, 1, Eigen::DontAlign> FtxPrime = fundamental_matrix.transpose() * x_cur;
    T JJt = ceres::sqrt(squ(Fx[0]) + squ(Fx[1]) + squ(FtxPrime[0]) + squ(FtxPrime[1]));

    return (x_cur.transpose() * Fx).eval()[0] / ceres::sqrt(JJt);
}
}

namespace torr_distance {
// as in "Invariant Fitting of Two View Geometry",Torr,Fitzgibbon
template <typename T>
T get_error(const Eigen::Matrix<T, 3, 3, Eigen::DontAlign>& fundamental_matrix,
            const Eigen::Matrix<T, 3, 1, Eigen::DontAlign>& x_cur,
            const Eigen::Matrix<T, 3, 1, Eigen::DontAlign>& x_old) {
    // calc sampson dist(Mulptiple view geom p.287)

    // nominator
    T xPrimetFx = x_cur.transpose() * fundamental_matrix * x_old;
    // denominator normalizes with the same effect as image_point normalization but also invariant
    // to equiform transformations
    T ftJf = fundamental_matrix.block(0, 0, 2, 2).squaredNorm();

    return xPrimetFx / ftJf;
}
}

namespace epipolar_distance_inverse_direction {

// first part of error function
template <typename T>
T get_error(const Eigen::Matrix<T, 3, 3, Eigen::DontAlign>& fundamental_matrix,
            const Eigen::Matrix<T, 3, 1, Eigen::DontAlign>& x_cur,
            const Eigen::Matrix<T, 3, 1, Eigen::DontAlign>& x_old) {
    // second part of symmetric epipolar error(Mulptiple view
    // geom(http://www.robots.ox.ac.uk/~az/tutorials/tutoriala.pdf))
    // to use it you have to add both epipolar_distance::get_error and this one as residual blocks
    // to the ceres problem
    Eigen::Matrix<T, 3, 1, Eigen::DontAlign> FtxPrime = fundamental_matrix.transpose() * x_cur;

    return (x_cur.transpose() * fundamental_matrix * x_old).eval()[0] /
           ceres::sqrt(squ(FtxPrime[0]) + squ(FtxPrime[1]));
}
}
namespace epipolar_distance {
template <typename T>
T get_error(const Eigen::Matrix<T, 3, 3, Eigen::DontAlign>& fundamental_matrix,
            const Eigen::Matrix<T, 3, 1, Eigen::DontAlign>& x_cur,
            const Eigen::Matrix<T, 3, 1, Eigen::DontAlign>& x_old) {
    // get epipolar line
    Eigen::Matrix<T, 3, 1, Eigen::DontAlign> Fx = fundamental_matrix * x_old;
    // normalize it
    return (x_cur.transpose() * Fx).eval()[0] / ceres::sqrt(squ(Fx[0]) + squ(Fx[1]));
    //    return (x_cur.transpose() * Fx).eval()[0]; // this error is used for 8 point since we
    //    equalize to zero -> no scaling necessary
}
}

namespace ray_to_epipolar_plane_linear {
template <typename T>
T get_error(const Eigen::Matrix<T, 3, 3, Eigen::DontAlign>& essential_matrix,
            const Eigen::Matrix<T, 3, 1, Eigen::DontAlign>& viewing_ray_cur,
            const Eigen::Matrix<T, 3, 1, Eigen::DontAlign>& viewing_ray_old) {
    // viewing_ray=v
    // epipolar constraint can be formulated as v_cur^T*essential_matrix*v_old=0. Since we use an
    // SVP model this means that essential_matrix*v_old=normal_epipolar_plane.
    // Therefore when v_cur.norm()==1. and essential_matrix*v_old==1. this error cos(a)
    // where "a" is the deviation angle of the epipolar plane+pi/2.
    // this error is approx. linear(since cos(x)|pi/2 ~=0-x+x^3/3!-...)
    Eigen::Matrix<T, 3, 1, Eigen::DontAlign> Ev = essential_matrix * viewing_ray_old;
    // normalize it
    Ev.normalize();

    auto sin_alpha = (viewing_ray_cur.transpose() / viewing_ray_cur.norm() * Ev).eval()[0];

    // scale the error, otherwise it will be too small
    return sin_alpha;
}
}

namespace ray_to_epipolar_plane_quadratic {
template <typename T>
T get_error(const Eigen::Matrix<T, 3, 3, Eigen::DontAlign>& essential_matrix,
            const Eigen::Matrix<T, 3, 1, Eigen::DontAlign>& viewing_ray_cur,
            const Eigen::Matrix<T, 3, 1, Eigen::DontAlign>& viewing_ray_old) {
    // viewing_ray=v
    // epipolar constraint can be formulated as v_cur^T*essential_matrix*v_old=0. Since we use an
    // SVP model this means that essential_matrix*v_old=normal_epipolar_plane.
    // Therefore when v_cur.norm()==1. and essential_matrix*v_old==1. this error equals sin(a) where
    // "a" is the deviation angle of the epipolar plane.

    // cosine is better conditioned around 0 (since cos(x)~= 1-x^2/2!+x^4/4!-... around zero)
    // therefore we take the projectino of the viewing ray onto the epipolar plane and calculate the
    // cosine of the angle between the plane and the ray

    // calculate normal of epipolar plane
    Eigen::Matrix<T, 3, 1, Eigen::DontAlign> Ev = essential_matrix * viewing_ray_old;
    // normalize it
    Ev.normalize();

    // calculate vector in epipolar plane
    Eigen::Matrix<T, 3, 1, Eigen::DontAlign> y = viewing_ray_cur - Ev.dot(viewing_ray_cur) * Ev;

    T cos_beta = viewing_ray_cur.dot(y) / y.norm() / viewing_ray_cur.norm();

    // scale the error, otherwise it will be too small
    return (static_cast<T>(1.) - cos_beta);
}
}
}
}
