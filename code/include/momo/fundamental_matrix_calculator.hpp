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

#include <memory>
#include <Eigen/Eigen>
#include "essential_matrix_calculator.hpp"

namespace momo {

/** @class FundamentalMatrixCalculator
 *  Wrap calibration aournd the essential matrix to get the fundamental matrix of a calibrated
camera, given a motion in a reference frame (f.e.
 * vehicle frame)
 * This is a one-way-class use it once than throw it away
 * When using Eigen with jet be carefull: from ceres/jet.h:
* // The infinitesimal part.
  //
  // Note the Eigen::DontAlign bit is needed here because this object
  // gets allocated on the stack and as part of other arrays and
  // structs. Forcing the right alignment there is the source of much
  // pain and suffering. Even if that works, passing Jets around to
  // functions by value has problems because the C++ ABI does not
  // guarantee alignment for function arguments.
  //
  // Setting the DontAlign bit prevents Eigen from using SSE for the
  // various operations on Jets. This is a small performance penalty
  // since the AutoDiff code will still expose much of the code as
  // statically sized loops to the compiler. But given the subtle
  // issues that arise due to alignment, especially when dealing with
  // multiple platforms, it seems to be a trade off worth making.
  Eigen::Matrix<T, N, 1, Eigen::DontAlign> v;
     */
template <typename T>
class FundamentalMatrixCalculator {
public:  // custom types
private: // attributes
    std::shared_ptr<const Camera> cam_;
    ceres_eigen_types::Matrix3<T> data;

public: // methods
    FundamentalMatrixCalculator(ceres_eigen_types::Affine3<T> motion_ref_frame,
                                std::shared_ptr<const Camera> cam)
            : cam_(cam) {
        this->data = Process(T(), motion_ref_frame);
    }

    ceres_eigen_types::Matrix3<T> get() const {
        return this->data;
    }

private:
    ///@brief calculates the fundamental matrix for the given camera in function of the motion in
    /// the reference frame
    ceres_eigen_types::Matrix3<T> Process(T, ceres_eigen_types::Affine3<T> motion_ref_frame) {
        using namespace ceres_eigen_types;

        EssentialMatrixCalculator<T> essential_matrix_calculator(
            motion_ref_frame, cam_); // instance for calculating the essential matrix

        // add calibration to get fundamental matrix
        Matrix3<T> fundamental_matrix =
            cam_->GetIntrinsicsInverse().transpose().cast<T>() * essential_matrix_calculator.get() *
            cam_->GetIntrinsicsInverse().cast<T>(); // multiple view geomtry p.257

        Assert(fundamental_matrix.allFinite(), "fundamental matrix has nan members");

        return fundamental_matrix;
    }

    void Assert(bool condition, std::string message) {
        if (!condition)
            throw std::runtime_error("in fundamental_matrix_calculator: " + message);
    }
};

} // end of ns momo
