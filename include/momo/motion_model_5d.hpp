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
#include <ceres/jet.h> // for ceres::sin
#include "motion_model_base.hpp"

namespace momo {

///@brief Reimplement Eigen::Angleaxis convertion since Eigen::AngleAxis cannot be used with
///Eigen::DontAlign which ceres' jet types need
template <typename T>
Eigen::Quaternion<T, Eigen::DontAlign> ToQuaternion(T angle, Eigen::Matrix<T, 3, 1,Eigen::DontAlign> axis) {
    // formula for quaternions
    // http://www.euclideanspace.com/maths/geometry/rotations/conversions/angleToQuaternion/
    auto static2 = static_cast<T>(2.);
    return Eigen::Quaternion<T, Eigen::DontAlign>(
        ceres::cos(angle / static2), axis[0] * ceres::sin(angle / static2),
        axis[1] * ceres::sin(angle / static2), axis[2] * ceres::sin(angle / static2));
}

/**
*  @class MotionModelOneWheel2d
*  @par
*
*  One wheel motion model: Movement on a plane along arclength l and yaw angle y
*  Implements the MotionModelBase interface
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
class MotionModel5d : public MotionModelBase<T> {
public:  // public classes/enums/types etc...
private: // attributes
    bool unset_;
    T arc_length_;
    T yaw_;

public: // public methods
    // default constructor
    MotionModel5d() : unset_(true) {
        ;
    }

    ///@brief constructor for setting input arguments arc length and yaw, for movement on a plane
    MotionModel5d(T arc_length, T yaw) : arc_length_(arc_length), yaw_(yaw), unset_(false) {
        ;
    }

    ///@brief wrapper for ceres sine for clear dependencies
    T sin(T a) const {
        return ceres::sin(a);
    }
    ///@brief wrapper for ceres cosine for clear dependencies
    T cos(T a) const {
        return ceres::cos(a);
    }

    ///@brief define custom assert with exception you can overload that with cassert if wanted
    void Assert(bool condition, std::string message) const {
        if (!condition)
            throw std::runtime_error("assertion happened in MotionModelOneWheel2d: " + message);
    }

    ///@brief calculate curvature without risk of dividing by zero
    T get_curvature() const {
        Assert(!IsUnset(), "model unset");
        T tmp = arc_length_;
        if (tmp == static_cast<T>(0.)) {
            tmp = static_cast<T>(0.001);
            std::cerr << "Attention arclength==0. -> curvature undefined" << std::endl;
        }

        return yaw_ / tmp;
    }

    T get_arc_length() const {
        return arc_length_;
    }

    T get_yaw() const {
        return yaw_;
    }

    bool IsUnset() const {
        return unset_;
    }

    ///@brief convert arc_length and yaw to 2d euclidean space (x,y)
    std::pair<T, T> ToEuclidean2d() const {
        Assert(!IsUnset(), "model unset");
        auto dx = arc_length_;
        auto dy = T(0.);

        auto curvature = get_curvature();
        if (curvature != T(0.)) {
            dx = sin(yaw_) / curvature;
            dy = (T(1.) - cos(yaw_)) / curvature;
        }
        return std::make_pair(dx, dy);
    }

    ///@brief interface for motion model
    Eigen::Transform<T, 3, Eigen::Affine, Eigen::DontAlign> ToPose() const override {
        using Affine3T = Eigen::Transform<T, 3, Eigen::Affine, Eigen::DontAlign>;
        using Vector3T = Eigen::Matrix<T, 3, 1, Eigen::DontAlign>;

        Assert(!IsUnset(), "model unset");
        T dy, dx;
        std::tie(dx, dy) = ToEuclidean2d();

        Affine3T out;
        Vector3T z_axis(static_cast<T>(0.), static_cast<T>(0.), static_cast<T>(1.));
        Eigen::Quaternion<T, Eigen::DontAlign> quat = ToQuaternion(yaw_, z_axis);
        out.linear() = quat.toRotationMatrix();
        out.translation() = Vector3T(dx, dy, static_cast<T>(0.));

        return out;
    }
};

} // end of ns
