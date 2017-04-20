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
#include <string>
#include <Eigen/Eigen>

namespace momo {

/** @class Camera
 *  Class for storing the camera props such as extrinsic and intrinsic calibration.
 *  Checks their validity.
 * */
class Camera {
public:
    // define pointer
    using Ptr = std::shared_ptr<Camera>;
    using ConstPtr = std::shared_ptr<const Camera>;

    Camera() {
        ;
    }

    ///@brief set id if you want to distinguish different cameras
    void set_id(std::string s) {
        id_ = s;
    }

    std::string get_id() const {
        return id_;
    }

    ///@brief constructor for non-pinhole models
    // extrinsics are defined as transformation from camera to reference frame
    Camera(Eigen::Affine3d pose_camera_to_reference)
            : pose_camera_to_reference_(pose_camera_to_reference) {
        pose_camera_to_reference_inverse_ = pose_camera_to_reference_.inverse();
    }

    ///@brief constructor for pinhole models
    // extrinsics are defined as transformation from camera to reference frame
    Camera(Eigen::Affine3d pose_camera_to_reference, Eigen::Matrix3d intrinsics)
            : pose_camera_to_reference_(pose_camera_to_reference), intrinsics_(intrinsics) {
        pose_camera_to_reference_inverse_ = pose_camera_to_reference_.inverse();
        intrinsics_inverse_ = intrinsics.inverse();
    }

    Eigen::Affine3d GetPoseCameraToReference() const {
        Assert(std::abs(pose_camera_to_reference_.rotation().determinant() - 1.) < 0.01,
               "rotation invalid");
        return pose_camera_to_reference_;
    }

    Eigen::Affine3d GetPoseReferenceToCamera() const {
        Assert(std::abs(pose_camera_to_reference_.rotation().determinant() - 1.) < 0.01,
               "rotation invalid");
        Assert(pose_camera_to_reference_inverse_.matrix().allFinite(), "divided by zero");
        return pose_camera_to_reference_inverse_;
    }

    Eigen::Matrix3d GetIntrinsics() const {
        Assert(intrinsics_(0, 0) > 1e-50, "intrinsics invalid");
        return intrinsics_;
    }

    Eigen::Matrix3d GetIntrinsicsInverse() const {
        Assert(intrinsics_(0, 0) > 1e-50, "intrsinsics invalid");
        Assert(intrinsics_inverse_.allFinite(), "divided by zero");
        return intrinsics_inverse_;
    }

private:
    ///@brief extrinsics_ is transformation from camera to reference frame (f.e. vehicle coordinate
    /// system)
    Eigen::Affine3d pose_camera_to_reference_, pose_camera_to_reference_inverse_;
    ///@brief intrinsics of the camera
    Eigen::Matrix3d intrinsics_, intrinsics_inverse_;
    ///@brief string to identify camera
    std::string id_ = "no id";

    ///@brief defines custom assert with exception. You can overload it with assert or custom
    /// exception if you want.
    void Assert(bool condition, std::string message) const {
        if (!condition)
            throw std::runtime_error("Assertion in Camera: " + message);
    }
};
}
