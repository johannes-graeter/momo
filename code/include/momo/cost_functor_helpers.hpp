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

#include <ceres/ceres.h>

#include "commons.hpp"
#include "cost_functors_ceres.hpp"

namespace momo {
namespace cost_functor_helpers {

// add residuals to ceres problem
// more information on the distances used can be found here http://www.robots.ox.ac.uk/~az/tutorials/tutoriala.pdf)) and in "Mulptiple View Geometry",Hartley and Zissermann, p.287)

void AddResidualSampsonDistance(ceres::Problem& problem,
                                const std::pair<Eigen::Vector3d, Eigen::Vector3d>& match_prev_cur,
                                double fix_l, std::shared_ptr<const Camera> camera) {
    ceres::CostFunction* cost_function = CostFunctorCeresRollPitchYaw::Create(
        match_prev_cur, fix_l, camera,
        commons::CostFunctionType::SampsonDistance);
    problem.AddResidualBlock(cost_function, NULL /** @todo add loss function option*/);
}

void AddResidualSymmetricEpipolarDistance(
    ceres::Problem& problem, const std::pair<Eigen::Vector3d, Eigen::Vector3d>& match_prev_cur,
    double fix_l, std::shared_ptr<const Camera> camera) {

    // add forward epipolar distance and backward epipolar distance as seperate residual blocks
    ceres::CostFunction* cost_function_forward = CostFunctorCeresRollPitchYaw::Create(
        match_prev_cur, fix_l, camera,
        commons::CostFunctionType::EpipolarDistanceForward);
    problem.AddResidualBlock(cost_function_forward, NULL /** @todo add loss function option*/);

    ceres::CostFunction* cost_function_backward = CostFunctorCeresRollPitchYaw::Create(
        match_prev_cur, fix_l, camera,
        commons::CostFunctionType::EpipolarDistanceBackward);
    problem.AddResidualBlock(cost_function_backward, NULL /** @todo add loss function option*/);
}

void AddResidualEpipolarDistance(
    ceres::Problem& problem, const std::pair<Eigen::Vector3d, Eigen::Vector3d>& match_prev_cur,
    double fix_l, std::shared_ptr<const Camera> camera) {

    // add forward epipoalr distacne and backward epipolar distance as seperate residual blocks
    ceres::CostFunction* cost_function_forward = CostFunctorCeresRollPitchYaw::Create(
        match_prev_cur, fix_l, camera,
        commons::CostFunctionType::EpipolarDistanceForward);
    problem.AddResidualBlock(cost_function_forward, NULL /** @todo add loss function option*/);
}
}
