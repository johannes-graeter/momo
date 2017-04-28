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

/**
*  @class MotionModelBase
*  Using google c++ style guide naming conventions
*
* interface for motion models that transfer optimization parameters to an affine3d
* set input parameters in  constructor
* templated because of Jet type for ceres optimization
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
class MotionModelBase {
public: // public classes/enums/types etc...
    using Affine3T = Eigen::Transform<T, 3, Eigen::Affine,Eigen::DontAlign>;
    using Vector3T = Eigen::Matrix<T, 3, 1,Eigen::DontAlign>;

public: // attributes
public: // public methods
    ///@brief default constructor
    MotionModelBase() {
        ;
    }

    ///@brief convert input to affine3
    virtual Affine3T ToPose() const = 0;
};
}
