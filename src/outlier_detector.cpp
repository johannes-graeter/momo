/*
 * Copyright 2017. All rights reserved.
 * Institute of Measurement and Control Systems
 * Karlsruhe Institute of Technology, Germany
 *
 * authors:
 *  Johannes Graeter (johannes.graeter@kit.edu)
 *  and others
 */
#include "outlier_detector.hpp"

namespace momo {

void OutlierDetectorBase::Assert(bool condition, std::string message) const {
    if (!condition)
        throw std::runtime_error("Assertion in OutlierDetectorBase: " + message);
}

OutlierDetectorThreshold::OutlierDetectorThreshold(double squared_threshold)
        : squared_threshold(squared_threshold) {
    ;
}

void OutlierDetectorThreshold::Assert(bool condition, std::string message) const {
    if (!condition)
        throw std::runtime_error("Assertion in OutlierDetectorThreshold: " + message);
}

OutlierDetectorBase::OutlierFlags OutlierDetectorThreshold::process(
    MotionOptimizer::MatchesConstPtr data) {
    OutlierFlags out;
    out.reserve(data->size());

    for (const auto& m : *data) {
        out.push_back(m.get_loss_squared() > squared_threshold);
    }
    return out;
}

void OutlierDetectorThreshold::set_squared_threshold(double a) {
    squared_threshold = a;
}

OutlierDetectorQuantile::OutlierDetectorQuantile(double quantile)
        : OutlierDetectorThreshold(-1.), quantile(quantile) {
    ;
}

void OutlierDetectorQuantile::Assert(bool condition, std::string message) const {
    if (!condition)
        throw std::runtime_error("Assertion in OutlierDetectorQuantile: " + message);
}

OutlierDetectorBase::OutlierFlags OutlierDetectorQuantile::process(
    MotionOptimizer::MatchesConstPtr data) {

    // sorting doesn't work if is zero return directly empty
    if(data->size()==0)
        return OutlierDetectorBase::OutlierFlags();

    // get quantile
    double quantile_thres;
    {
        std::vector<double> residuals; // how much does it take to copy that?
        residuals.reserve(data->size());

        for (const auto& el : *data) {
            residuals.push_back(el.get_loss_squared());
        }

        // sort so that neth element is in correct place
        int index = static_cast<int>(quantile * (residuals.size() - 1));
        std::nth_element(residuals.begin(), std::next(residuals.begin(), index), residuals.end());
        quantile_thres = residuals[index];
    }

    // since quantile thresholding is a special kind of fix thresholding call parent function
    this->set_squared_threshold(quantile_thres);
    return OutlierDetectorThreshold::process(data);
}
}
