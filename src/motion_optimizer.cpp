/*
 * Copyright 2017. All rights reserved.
 * Institute of Measurement and Control Systems
 * Karlsruhe Institute of Technology, Germany
 *
 * authors:
 *  Johannes Graeter (johannes.graeter@kit.edu)
 *  and others
 */
#include <cassert>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#include <ceres/ceres.h>
#include "cost_functors_ceres.hpp"

#include "motion_optimizer.hpp"

namespace momo {

// estimation problem doesn't move from [0.,0.,0.] because R*[t]_x is all zero for any t if
// R==Identity, set it to a different value
// MotionOptimizer::MotionOptimizer() : last_estimate_{0.05} {
MotionOptimizer::MotionOptimizer() : last_estimate_{0.05, 0.05, 0.05} {
    ;
}

void MotionOptimizer::SetSolverOptions(MotionOptimizer::Options o) {
    solver_options_ = o;
}

MotionOptimizer::Options MotionOptimizer::GetSolverOptions() {
    return solver_options_;
}

void MotionOptimizer::ClearData() {
    data_.clear();
}

void MotionOptimizer::AddData(MotionOptimizer::MatchesPtr matches, Camera::ConstPtr camera,
                              commons::CostFunctionType type) {
    data_.push_back(std::make_tuple(matches, camera, type));
}

bool MotionOptimizer::IsReadyToOptimize() {
    bool is_ready = true;
    is_ready = is_ready && (data_.size() > 0);

    return is_ready;
}

// void MotionOptimizer::SetStartValues(std::vector<double> start_roll_pitch_yaw) {
//    last_estimate_ = start_roll_pitch_yaw;
//}

std::string MotionOptimizer::GetReport() {
    return report_;
}

// evaluate problem ones more and store residuals in data, use piotners for input args
void store_residuals(const std::vector<std::vector<ceres::ResidualBlockId>>& residual_block_ids,
                     ceres::Problem* problem, MotionOptimizer::StorageDataType* data) {
    // store residuals per match for in data_
    auto i1 = residual_block_ids.cbegin();
    auto i2 = data->begin();
    if (residual_block_ids.size() != data->size()) {
        throw std::runtime_error("residual_block_ids.size()=" +
                                 std::to_string(residual_block_ids.size()) + " != data.size()=" +
                                 std::to_string(data->size()));
    }

    // iterate over each camera and get the corresponding residuals
    // store them in data
    for (; i1 != residual_block_ids.cend() && i2 != data->end(); ++i1, ++i2) {
        // apparently ceres gives back other residuals when there is ids are empty
        // so this check is necessary
        if (i1->size() > 0) {
            // evaluate once with and once without loss function and get difference
            ceres::Problem::EvaluateOptions eval_options;
            eval_options.residual_blocks = *i1;
            // evaluate without loss function
            std::vector<double> residuals_without_loss;
            {
                eval_options.apply_loss_function = false;
                double cost = -1.;
                problem->Evaluate(eval_options, &cost, &residuals_without_loss, nullptr, nullptr);
            }
            std::vector<double> residuals_with_loss;
            {
                eval_options.apply_loss_function = true;
                double cost = -1.;
                problem->Evaluate(eval_options, &cost, &residuals_with_loss, nullptr, nullptr);
            }

            // get difference
            std::vector<double> residuals;
            residuals.reserve(residuals_without_loss.size());
            std::transform(residuals_without_loss.cbegin(), residuals_without_loss.cend(),
                           residuals_with_loss.cbegin(), std::back_inserter(residuals),
                           std::minus<double>());

            // assign them to matches in data
            // we have to be carefull that residuals are stored in same order as matches in data
            auto residuals_iter = residuals.cbegin();
            auto& matches = *std::get<0>(*i2); // for clarity
            auto matches_iter = matches.begin();
            std::cout << "residuals.size=" << residuals.size()
                      << " matches.size()=" << matches.size() << std::endl;
            if (residuals.size() != matches.size()) {
                throw std::runtime_error("residuals.size()=" + std::to_string(residuals.size()) +
                                         " != matches.size()=" + std::to_string(matches.size()));
            }
            for (; residuals_iter != residuals.cend() && matches_iter != matches.end();
                 ++residuals_iter, ++matches_iter) {
                matches_iter->set_loss(*residuals_iter);
            }
        }
    }
}

std::string MotionOptimizer::Solve(std::vector<double>& x, double fix_l, bool verbose) {

    //    auto start_prep_solve = std::chrono::high_resolution_clock::now();

    ceres::Problem problem;

    // store residual block ids, so that we can evaluate the residuals once more at the end and add
    // them to matches
    std::vector<std::vector<ceres::ResidualBlockId>> residual_block_ids;
    residual_block_ids.reserve(data_.size());

    double scale = 1e+10;

    // add residuals
    for (const auto& el : data_) {
        MatchesConstPtr matches;
        Camera::ConstPtr cam;
        commons::CostFunctionType cost_type;

        std::tie(matches, cam, cost_type) = el;

        if (verbose)
            std::cout << "adding " << matches->size()
                      << " matches to the problem.\n fix_l=" << fix_l
                      << "\npose camera to reference=\n"
                      << cam->GetPoseCameraToReference().matrix() << "\ncost_type="
                      << static_cast<std::underlying_type<commons::CostFunctionType>::type>(
                             cost_type)
                      << std::endl;

        // add residual block with loss function, lossfunctions are described in
        // https://github.com/kashif/ceres-solver/blob/master/include/ceres/loss_function.h
        //        ceres::ScaledLoss* loss_function(new ceres::CauchyLoss(0.5), 10000.,
        //        ceres::TAKE_OWNERSHIP);

        // store residual block ids
        std::vector<ceres::ResidualBlockId> cur_ids;
        cur_ids.reserve(matches->size());

        for (const auto& match : *matches) {
            auto cur_id = problem.AddResidualBlock(
                CostFunctorCeresRollPitchYaw::Create(match.MakePair(), fix_l, cam, cost_type),
                new ceres::ScaledLoss(new ceres::CauchyLoss(solver_options_.loss_function_width),
                                      scale, ceres::TAKE_OWNERSHIP),
                x.data());
            cur_ids.push_back(std::move(cur_id));
        }

        residual_block_ids.push_back(std::move(cur_ids));
    }

    //    if (x.size() == 1) {
    //        problem.SetParameterLowerBound(x.data(), 0, -1.);
    //        problem.SetParameterUpperBound(x.data(), 0, 1.);
    //    }

    //    auto end_prep_solve = std::chrono::high_resolution_clock::now();
    //    int64_t duration_prep_solve =
    //        std::chrono::duration_cast<std::chrono::milliseconds>(end_prep_solve -
    //        start_prep_solve)
    //            .count();
    //    if (verbose)
    //        std::cout << "Duration prep of MotionOptimizer solve: " << duration_prep_solve << "
    //        ms"
    //                  << std::endl;

    // Run the solver!
    ceres::Solver::Options options;
    options.linear_solver_type = ceres::DENSE_QR;
    options.max_num_iterations = solver_options_.max_number_iterations;
    options.max_solver_time_in_seconds = solver_options_.max_time_solving_sec;
    options.function_tolerance = solver_options_.function_tolerance;
    options.gradient_tolerance = solver_options_.gradient_tolerance;
    options.parameter_tolerance = solver_options_.parameter_tolerance;

    if (verbose)
        options.minimizer_progress_to_stdout = true;

    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);

    store_residuals(residual_block_ids, &problem, &data_);

    relative_cost_change_ = (summary.initial_cost - summary.final_cost) / summary.initial_cost;

    return summary.FullReport();
}

Eigen::Affine3d MotionOptimizer::Solve(double arc_length, bool verbose) {

    // get estimate, store it in last_estimate_roll_pitch_yaw
    if (IsReadyToOptimize()) {
        report_ = Solve(last_estimate_, arc_length, verbose);
    } else {
        report_ = "not ready to optimize";
    }

    //    if (verbose)
    //        std::cout << "params after optimization estimate: roll=" << last_estimate_[0]
    //                  << " pitch=" << last_estimate_[1] << " yaw=" << last_estimate_[2] <<
    //                  std::endl;

    // apply motion model
    MotionModelOneWheel3d<double> motion_model(arc_length, last_estimate_[0], last_estimate_[1],
                                               last_estimate_[2]);

    //    // for high yaw problem can be stuck on a plane surface in the error landscape -> reset to
    //    0.05 yaw
    if (relative_cost_change_ < solver_options_.relative_cost_change_to_reset_prior) {
        last_estimate_ = std::vector<double>{0.05, 0.05, 0.05};
        report_ = report_ + "; reset roll, pitch, yaw to 0.05 rad";
    }

    return motion_model.ToPose();
}
}
