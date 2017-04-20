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
#include "ceres_eigen_types.hpp"
#include "commons.hpp"

namespace momo {
/** @class EssentialMatrixCalculator
 *  Get the essential matrix of the epipoar geometry, given a motion in a reference frame (f.e.
 * vehicle frame)
 * This is a one-way-class use it once then throw it away
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
class EssentialMatrixCalculator {
public:  // custom types
private: // attributes
    std::shared_ptr<const Camera> cam_;
    ceres_eigen_types::Matrix3<T> data;

public: // methods
    EssentialMatrixCalculator(ceres_eigen_types::Affine3<T> motion_ref_frame,
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

        Assert(motion_ref_frame.matrix().allFinite(), "motion_ref_frame has nan");

        // transfer motion
        Affine3<T> motion_cam = CalculateCameraMotion(motion_ref_frame);

        Assert(motion_cam.matrix().allFinite(), "motion_cam has nan");

        // calculate essential matrix
        Matrix3<T> essential_matrix =
            motion_cam.rotation().inverse() * commons::GetCrossProdMatrix(motion_cam.translation().eval());


        Assert(essential_matrix.allFinite(), "essential matrix has nan members");

        return essential_matrix;
    }

    ///@brief transfers the  motion of the reference frame into the camera frame
    ceres_eigen_types::Affine3<T> CalculateCameraMotion(
        ceres_eigen_types::Affine3<T> reference_motion) {
        // this is the rear camera motion 3d, not the motion of points in the camera frame
        ceres_eigen_types::Affine3<T> camera_motion = cam_->GetPoseCameraToReference().cast<T>() *
                                                      reference_motion *
                                                      cam_->GetPoseReferenceToCamera().cast<T>();
        // scale of translation is a scalar factor in the optimization -> better convergence if
        // normalized?
        // for very small motion this gets messy
        //        camera_motion.translation().normalize();
        return camera_motion;
    }

    void Assert(bool condition, std::string message) {
        if (!condition)
            throw std::runtime_error("in essential_matrix_calculator: " + message);
    }
};
}
