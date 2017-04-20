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
#include "motion_optimizer.hpp"

namespace momo {

/**
*  @class OutlierDetectorBase
*
*  base class for outlier detectors that take the loss(difference of resiudals with and without loss
* function) and decide if corresponding match is an
* outlier or not
*
*  Uses google c++ style guide naming conventions
*/
class OutlierDetectorBase {
public: // public classes/enums/types etc...
    using OutlierFlags =
        std::vector<bool>; ///< if true corresponding match is outlier, order matters!

private: // attributes
public:  // public methods
    // default constructor
    OutlierDetectorBase() {
        ;
    }

    virtual OutlierFlags process(MotionOptimizer::MatchesConstPtr data) = 0;

    ///@brief defines custom assert with exception. You can overload it with assert if you want
    virtual void Assert(bool condition, std::string message) const;
}; // end of OutlierDetectorBase

/**
*  @class OutlierDetectorThreshold
*
*  classic approach: use fix threshold to decide wether its an inlier or outlier
*
*  Uses google c++ style guide naming conventions
*/
class OutlierDetectorThreshold : public OutlierDetectorBase {
public:                       // public classes/enums/types etc...
private:                      // attributes
    double squared_threshold; ///< fix threshold
public:                       // public methods
    // default constructor
    OutlierDetectorThreshold(double squared_threshold = -1.);

    ///@brief defines custom assert with exception. You can overload it with assert if you want
    virtual void Assert(bool condition, std::string message) const override;

    OutlierFlags process(MotionOptimizer::MatchesConstPtr data) override;

    void set_squared_threshold(double a);
}; // end of OutlierDetectorThreshold

/**
*  @class OutlierDetectorQuantile
*
*  threshold is relative to residual value: if quantile=0.9, the upper 10% of residuals will be
* rejected
*
*  Uses google c++ style guide naming conventions
*/
class OutlierDetectorQuantile : public OutlierDetectorThreshold {
public:              // public classes/enums/types etc...
private:             // attributes
    double quantile; ///< quantile, if =0.9 upper 10% of measurement will be rejected
public:              // public methods
    // default constructor
    OutlierDetectorQuantile(double quantile);

    ///@brief defines custom assert with exception. You can overload it with assert if you want
    virtual void Assert(bool condition, std::string message) const override;

    OutlierFlags process(MotionOptimizer::MatchesConstPtr data) override;
};

// end of OutlierDetectorQuantile
}
