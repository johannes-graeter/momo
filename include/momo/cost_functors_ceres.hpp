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

//#include <cassert>
//#include <cmath>
//#include <iomanip>
//#include <iostream>
//#include <map>
//#include <sstream>
//#include <string>
//#include <tuple>
//#include <vector>

#include <memory>

#include <Eigen/Eigen>
#include <ceres/ceres.h>

///@todo make camera calibration container
#include "camera.hpp"

#include "commons.hpp"
#include "cost_functions.hpp"
#include "essential_matrix_calculator.hpp"
#include "fundamental_matrix_calculator.hpp"
#include "motion_model_one_wheel.hpp"

namespace momo {

/**
*  @class CostFunctorCeresRollPitchYaw
*
*  Cost functor for one single flow vector for usage in google's ceres solver.
*  Different error functions can be defined by a flag in the constructor, set with void
set_cost_function_type(CostFunctionType type)
*  Use in ceres:
*  ceres::CostFunction* cost_function = new
* ceres::AutoDiffCostFunction<momo::CostFunctorCeresRollPitchYaw, 1, 3>(new
* CostFunctorCeresRollPitchYaw(std::make_pair(match_previous,match_current), arc_length, camera);
*
* or use the ::Create function
*
* for adding residuals refer to cost_functor_helpers.hpp
*
* You find a basic example in http://ceres-solver.org/nnls_tutorial.html
*
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
class CostFunctorCeresRollPitchYaw {
public:                          // public classes/enums/types etc...
    using Vec = Eigen::Vector3d; // matches get stored as vectors in 3d in most cases Vec a[2]==1.

    struct ExceptionZeroArcLength : public std::exception {
        virtual const char* what() const throw() {
            return "in CostFunctorCeresRollPitchYaw: fixed arclength is <=0 -> cost functor is "
                   "undefined";
        }
    };

public: // attributes
    ///@brief here we do not optimize the arclength so the arclength has to be determined from
    /// outside
    /// @todo add variant with arc length as a variable
    double arc_length_fixed_;
    ///@brief the two matches we want to have the error from
    Vec match_previous_, match_current_;
    ///@brief storage for calibration and so on
    std::shared_ptr<const Camera> cam_;
    ///@brief tell the compiler which error function to use
    commons::CostFunctionType cost_function_type_;

public: // public methods
    // default constructor
    CostFunctorCeresRollPitchYaw(const std::pair<Vec, Vec>& match_prev_cur, double fix_l,
                                 std::shared_ptr<const Camera> camera,
                                 commons::CostFunctionType type = commons::CostFunctionType::SampsonDistance)
            : arc_length_fixed_(fix_l), cam_(camera), match_previous_(match_prev_cur.first),
              match_current_(match_prev_cur.second), cost_function_type_(type) {
        if (arc_length_fixed_ < 0.001) {
            throw ExceptionZeroArcLength();
        }
    }

    void set_cost_function_type(commons::CostFunctionType type) {
        cost_function_type_ = type;
    }

    ///@brief the matches and the fundamental matrix evaluate to this error (f.e. sampson
    /// distance,(symmetric) epipolar distance,...)

    // this operator gets called by ceres in each iteration step several times (evaluation and
    // gradient calculation) parameters are x=[roll,pitch,yaw] in vreference frame coordinates
    template <typename T>
    bool operator()(const T* const x, T* residual) const {
        using namespace commons;

        using Affine3T = ceres_eigen_types::Affine3<T>;
        using Matrix3T = ceres_eigen_types::Matrix3<T>;
        using Vector3T = ceres_eigen_types::Vector3<T>;

        // transform input as x=[roll, pitch, yaw] and fix arclength to a motion of the reference
        // frame

        MotionModelOneWheel3d<T> motion_model(T(arc_length_fixed_), x[0], x[1], x[2]);
        Affine3T motion_in_reference_frame = motion_model.ToPose();
        Assert(motion_in_reference_frame.matrix().allFinite(), "non finite motion");

        // get the fundamental matrix or the essential, depending on the omaine of the cost function
        // with one way class
        Matrix3T f;
        if (cost_function_type_ == CostFunctionType::RayToEpipolarPlaneLinear ||
            cost_function_type_ == CostFunctionType::RayToEpipolarPlaneQuadratic) {
            EssentialMatrixCalculator<T> f_calc(motion_in_reference_frame, cam_);
            f = f_calc.get();
        } else {
            FundamentalMatrixCalculator<T> f_calc(motion_in_reference_frame, cam_);
            f = f_calc.get();
        }

        Assert(f.matrix().allFinite(), "non finite fundamental matrix");

        // cast match to T (for jet type)
        Vector3T match_cur_t = match_current_.cast<T>();
        Vector3T match_prev_t = match_previous_.cast<T>();

        // calculate with cost function residual
        ///@todo do with pointer to cost functor -> curiously recurring tempalte pattern
        switch (cost_function_type_) {
        case CostFunctionType::SampsonDistance:
            residual[0] = cost_functions::sampson_distance::get_error(f, match_cur_t, match_prev_t);
            break;
        case CostFunctionType::EpipolarDistanceForward:
            residual[0] =
                cost_functions::epipolar_distance::get_error(f, match_cur_t, match_prev_t);
            break;
        case CostFunctionType::EpipolarDistanceBackward:
            residual[0] = cost_functions::epipolar_distance_inverse_direction::get_error(
                f, match_cur_t, match_prev_t);
            break;
        case CostFunctionType::RayToEpipolarPlaneLinear:
            residual[0] =
                cost_functions::ray_to_epipolar_plane_linear::get_error(f, match_cur_t, match_prev_t);
            break;
        case CostFunctionType::RayToEpipolarPlaneQuadratic:
            residual[0] =
                cost_functions::ray_to_epipolar_plane_quadratic::get_error(f, match_cur_t, match_prev_t);
            break;
        default:
            throw std::runtime_error("in CostFunctorCeresRollPitchYaw: undefined cost function "
                                     "type, must be one of {SampsonDistance, "
                                     "EpipolarDistanceForward, EpipolarDistanceBackward}");
            break;
        }

//        std::cout<<residual[0]<<std::endl;

        // set output and return
        return true;
    }

    // Factory to hide the construction of the CostFunction object from
    // the client code. Needs dependency from ceres.
    static ceres::CostFunction* Create(const std::pair<Vec, Vec>& match_prev_cur, double fix_l,
                                       Camera::ConstPtr camera,
                                       commons::CostFunctionType type = commons::CostFunctionType::SampsonDistance) {
        return (new ceres::AutoDiffCostFunction<CostFunctorCeresRollPitchYaw, 1, 3>(
            new CostFunctorCeresRollPitchYaw(match_prev_cur, fix_l, camera, type)));
    }

private:
    ///@brief defines custom assert with exception. You can overload it with assert or custom
    /// exception if you want.
    void Assert(bool condition, std::string message) const {
        if (!condition)
            throw std::runtime_error("Assertion in CostFunctorsCeres: " + message);
    }
};

/**
*  @class CostFunctorCeres5d
*
*  Cost functor for one single flow vector for usage in google's ceres solver.
*  Different error functions can be defined by a flag in the constructor, set with void
set_cost_function_type(CostFunctionType type)
*  Use in ceres:
*  ceres::CostFunction* cost_function = new
* ceres::AutoDiffCostFunction<momo::CostFunctorCeresRollPitchYaw, 1, 3>(new
* CostFunctorCeresRollPitchYaw(std::make_pair(match_previous,match_current), arc_length, camera);
*
* or use the ::Create function
*
* for adding residuals refer to cost_functor_helpers.hpp
*
* You find a basic example in http://ceres-solver.org/nnls_tutorial.html
*
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
/*
class CostFunctorCeres5d {
public:                          // public classes/enums/types etc...
    using Vec = Eigen::Vector3d; // matches get stored as vectors in 3d in most cases Vec a[2]==1.

    struct ExceptionZeroTranslation : public std::exception {
        virtual const char* what() const throw() {
            return "in CostFunctorCeres5d: translation is <=0 -> cost functor is "
                   "undefined";
        }
    };

public: // attributes
    ///@brief here we do not optimize the arclength so the arclength has to be determined from
    /// outside
    /// @todo add variant with arc length as a variable
    double translation_fixed_;
    ///@brief the two matches we want to have the error from
    Vec match_previous_, match_current_;
    ///@brief storage for calibration and so on
    std::shared_ptr<const Camera> cam_;
    ///@brief tell the compiler which error function to use
    commons::CostFunctionType cost_function_type_;

public: // public methods
    // default constructor
    CostFunctorCeresRollPitchYaw(const std::pair<Vec, Vec>& match_prev_cur, double fix_translation,
                                 std::shared_ptr<const Camera> camera,
                                 commons::CostFunctionType type = commons::CostFunctionType::SampsonDistance)
            : translation_fixed_(fix_translation), cam_(camera), match_previous_(match_prev_cur.first),
              match_current_(match_prev_cur.second), cost_function_type_(type) {
        if (translation_fixed_ < 0.001) {
            throw ExceptionZeroTranslation();
        }
    }

    void set_cost_function_type(commons::CostFunctionType type) {
        cost_function_type_ = type;
    }

    ///@brief the matches and the fundamental matrix evaluate to this error (f.e. sampson
    /// distance,(symmetric) epipolar distance,...)

    // this operator gets called by ceres in each iteration step several times (evaluation and
    // gradient calculation) parameters are x=[roll,pitch,yaw] in vreference frame coordinates
    template <typename T>
    bool operator()(const T* const x, T* residual) const {
        using namespace commons;

        using Affine3T = ceres_eigen_types::Affine3<T>;
        using Matrix3T = ceres_eigen_types::Matrix3<T>;
        using Vector3T = ceres_eigen_types::Vector3<T>;

        // transform input as x=[roll, pitch, yaw] and fix arclength to a motion of the reference
        // frame

        MotionModel5d<T> motion_model(T(translation_fixed_), x[0], x[1], x[2], x[3], x[4]);
        Affine3T motion_in_reference_frame = motion_model.ToPose();
        Assert(motion_in_reference_frame.matrix().allFinite(), "non finite motion");

        // get the fundamental matrix or the essential, depending on the omaine of the cost function
        // with one way class
        Matrix3T f;
        if (cost_function_type_ == CostFunctionType::RayToEpipolarPlaneLinear ||
            cost_function_type_ == CostFunctionType::RayToEpipolarPlaneQuadratic) {
            EssentialMatrixCalculator<T> f_calc(motion_in_reference_frame, cam_);
            f = f_calc.get();
        } else {
            FundamentalMatrixCalculator<T> f_calc(motion_in_reference_frame, cam_);
            f = f_calc.get();
        }

        Assert(f.matrix().allFinite(), "non finite fundamental matrix");

        // cast match to T (for jet type)
        Vector3T match_cur_t = match_current_.cast<T>();
        Vector3T match_prev_t = match_previous_.cast<T>();

        // calculate with cost function residual
        ///@todo do with pointer to cost functor -> curiously recurring tempalte pattern
        switch (cost_function_type_) {
        case CostFunctionType::SampsonDistance:
            residual[0] = cost_functions::sampson_distance::get_error(f, match_cur_t, match_prev_t);
            break;
        case CostFunctionType::EpipolarDistanceForward:
            residual[0] =
                cost_functions::epipolar_distance::get_error(f, match_cur_t, match_prev_t);
            break;
        case CostFunctionType::EpipolarDistanceBackward:
            residual[0] = cost_functions::epipolar_distance_inverse_direction::get_error(
                f, match_cur_t, match_prev_t);
            break;
        case CostFunctionType::RayToEpipolarPlaneLinear:
            residual[0] =
                cost_functions::ray_to_epipolar_plane_linear::get_error(f, match_cur_t, match_prev_t);
            break;
        case CostFunctionType::RayToEpipolarPlaneQuadratic:
            residual[0] =
                cost_functions::ray_to_epipolar_plane_quadratic::get_error(f, match_cur_t, match_prev_t);
            break;
        default:
            throw std::runtime_error("in CostFunctorCeresRollPitchYaw: undefined cost function "
                                     "type, must be one of {SampsonDistance, "
                                     "EpipolarDistanceForward, EpipolarDistanceBackward}");
            break;
        }

        // set output and return
        return true;
    }

    // Factory to hide the construction of the CostFunction object from
    // the client code. Needs dependency from ceres.
    static ceres::CostFunction* Create(const std::pair<Vec, Vec>& match_prev_cur, double fix_t,
                                       Camera::ConstPtr camera,
                                       commons::CostFunctionType type = commons::CostFunctionType::SampsonDistance) {
        return (new ceres::AutoDiffCostFunction<CostFunctorCeres5d, 1, 5>(
            new CostFunctorCeres5d(match_prev_cur, fix_t, camera, type)));
    }

private:
    ///@brief defines custom assert with exception. You can overload it with assert or custom
    /// exception if you want.
    void Assert(bool condition, std::string message) const {
        if (!condition)
            throw std::runtime_error("Assertion in CostFunctorsCeres: " + message);
    }
};
*/

///@todo make this more modular cause only motion model is different

/**
*  @class CostFunctorCeresYaw
*
*  Cost functor for one single flow vector for usage in google's ceres solver.
*  Different error functions can be defined by a flag in the constructor, set with void
set_cost_function_type(CostFunctionType type)
*  Use in ceres:
*  ceres::CostFunction* cost_function = new
* ceres::AutoDiffCostFunction<momo::CostFunctorCeresRollPitchYaw, 1, 3>(new
* CostFunctorCeresRollPitchYaw(std::make_pair(match_previous,match_current), arc_length, camera);
*
* or use the ::Create function
*
* for adding residuals refer to cost_functor_helpers.hpp
*
* You find a basic example in http://ceres-solver.org/nnls_tutorial.html
*
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
class CostFunctorCeresYaw {
public:                          // public classes/enums/types etc...
    using Vec = Eigen::Vector3d; // matches get stored as vectors in 3d in most cases Vec a[2]==1.

    struct ExceptionZeroArcLength : public std::exception {
        virtual const char* what() const throw() {
            return "in CostFunctorCeresRollPitchYaw: fixed arclength is <=0 -> cost functor is "
                   "undefined";
        }
    };

public: // attributes
    ///@brief here we do not optimize the arclength so the arclength has to be determined from
    /// outside
    /// @todo add variant with arc length as a variable
    double arc_length_fixed_;
    ///@brief the two matches we want to have the error from
    Vec match_previous_, match_current_;
    ///@brief storage for calibration and so on
    std::shared_ptr<const Camera> cam_;
    ///@brief tell the compiler which error function to use
    commons::CostFunctionType cost_function_type_;

public: // public methods
    // default constructor
    CostFunctorCeresYaw(const std::pair<Vec, Vec>& match_prev_cur, double fix_l,
                                 std::shared_ptr<const Camera> camera,
                                 commons::CostFunctionType type = commons::CostFunctionType::SampsonDistance)
            : arc_length_fixed_(fix_l), cam_(camera), match_previous_(match_prev_cur.first),
              match_current_(match_prev_cur.second), cost_function_type_(type) {
        if (arc_length_fixed_ < 0.001) {
            throw ExceptionZeroArcLength();
        }
    }

    void set_cost_function_type(commons::CostFunctionType type) {
        cost_function_type_ = type;
    }

    ///@brief the matches and the fundamental matrix evaluate to this error (f.e. sampson
    /// distance,(symmetric) epipolar distance,...)

    // this operator gets called by ceres in each iteration step several times (evaluation and
    // gradient calculation) parameters are x=[roll,pitch,yaw] in vreference frame coordinates
    template <typename T>
    bool operator()(const T* const x, T* residual) const {
        using namespace commons;

        using Affine3T = ceres_eigen_types::Affine3<T>;
        using Matrix3T = ceres_eigen_types::Matrix3<T>;
        using Vector3T = ceres_eigen_types::Vector3<T>;

        // transform input as x=[roll, pitch, yaw] and fix arclength to a motion of the reference
        // frame

        MotionModelOneWheel2d<T> motion_model(T(arc_length_fixed_), x[0]);
        Affine3T motion_in_reference_frame = motion_model.ToPose();
        Assert(motion_in_reference_frame.matrix().allFinite(), "non finite motion");

        // get the fundamental matrix or the essential, depending on the omaine of the cost function
        // with one way class
        Matrix3T f;
        if (cost_function_type_ == CostFunctionType::RayToEpipolarPlaneLinear ||
            cost_function_type_ == CostFunctionType::RayToEpipolarPlaneQuadratic) {
            EssentialMatrixCalculator<T> f_calc(motion_in_reference_frame, cam_);
            f = f_calc.get();
        } else {
            FundamentalMatrixCalculator<T> f_calc(motion_in_reference_frame, cam_);
            f = f_calc.get();
        }

        Assert(f.matrix().allFinite(), "non finite fundamental matrix");

        // cast match to T (for jet type)
        Vector3T match_cur_t = match_current_.cast<T>();
        Vector3T match_prev_t = match_previous_.cast<T>();

        // calculate with cost function residual
        ///@todo do with pointer to cost functor -> curiously recurring tempalte pattern
        switch (cost_function_type_) {
        case CostFunctionType::SampsonDistance:
            residual[0] = cost_functions::sampson_distance::get_error(f, match_cur_t, match_prev_t);
            break;
        case CostFunctionType::EpipolarDistanceForward:
            residual[0] =
                cost_functions::epipolar_distance::get_error(f, match_cur_t, match_prev_t);
            break;
        case CostFunctionType::EpipolarDistanceBackward:
            residual[0] = cost_functions::epipolar_distance_inverse_direction::get_error(
                f, match_cur_t, match_prev_t);
            break;
        case CostFunctionType::RayToEpipolarPlaneLinear:
            residual[0] =
                cost_functions::ray_to_epipolar_plane_linear::get_error(f, match_cur_t, match_prev_t);
            break;
        case CostFunctionType::RayToEpipolarPlaneQuadratic:
            residual[0] =
                cost_functions::ray_to_epipolar_plane_quadratic::get_error(f, match_cur_t, match_prev_t);
            break;
        default:
            throw std::runtime_error("in CostFunctorCeresRollPitchYaw: undefined cost function "
                                     "type, must be one of {SampsonDistance, "
                                     "EpipolarDistanceForward, EpipolarDistanceBackward}");
            break;
        }

        // set output and return
        return true;
    }

    // Factory to hide the construction of the CostFunction object from
    // the client code. Needs dependency from ceres.
    static ceres::CostFunction* Create(const std::pair<Vec, Vec>& match_prev_cur, double fix_l,
                                       Camera::ConstPtr camera,
                                       commons::CostFunctionType type = commons::CostFunctionType::SampsonDistance) {
        return (new ceres::AutoDiffCostFunction<CostFunctorCeresYaw, 1, 1>(
            new CostFunctorCeresYaw(match_prev_cur, fix_l, camera, type)));
    }

private:
    ///@brief defines custom assert with exception. You can overload it with assert or custom
    /// exception if you want.
    void Assert(bool condition, std::string message) const {
        if (!condition)
            throw std::runtime_error("Assertion in CostFunctorsCeres: " + message);
    }
};

}
