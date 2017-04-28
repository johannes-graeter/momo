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

#include "camera.hpp"
//#include <ceres/ceres.h>
#include "commons.hpp" //for cost function type

namespace momo {
/**
*  @class MotionOptimizer
*  @par Wrapper for ceres optimization problem
*   minimizes one of the cost functors defined in cost_functors_ceres.hpp
*
*/
class MotionOptimizer {
public:
    /// ///////////////////////////////////////////
    ///@brief options for ceres problem, is a subset of ceres options. Don't want to have dep of
    /// ceres in Header
    struct Options {

        Options() {
            max_number_iterations = 50;
            max_time_solving_sec = 1e6;
            function_tolerance = 1e-6;
            gradient_tolerance = 1e-10;
            parameter_tolerance = 1e-8;
            relative_cost_change_to_reset_prior = 1e-10;
            loss_function_width = 0.1;
        }

        int max_number_iterations;
        double max_time_solving_sec;
        double function_tolerance;
        double gradient_tolerance;
        double parameter_tolerance;
        double relative_cost_change_to_reset_prior;
        double loss_function_width;

        friend std::ostream& operator<<(std::ostream& stream, const Options& data) {
            stream << "max_number_iterations=" << data.max_number_iterations << std::endl;
            stream << "max_time_solving_sec=" << data.max_time_solving_sec << std::endl;
            stream << "function_tolerance=" << data.function_tolerance << std::endl;
            stream << "gradient_tolerance=" << data.gradient_tolerance << std::endl;
            stream << "parameter_tolerance=" << data.parameter_tolerance << std::endl;
            stream << "relative_cost_change_to_reset_prior="
                   << data.relative_cost_change_to_reset_prior << std::endl;
            stream << "loss_function_width=" << data.loss_function_width << std::endl;

            return stream;
        }
    };

    ///@brief u,v,1 or viewing ray(for geometric errors it doesn't matter how viewing ray is normed,
    /// will be normalized in cost functor
    using Ray = Eigen::Vector3d;

    /// ////////////////////////////////////////////
    /// @brief Match - Container for matches
    struct Match {
        Match() {
            ;
        }
        Match(Ray p_prev, Ray p_cur) : previous(p_prev), current(p_cur) {
            ;
        }
        Ray previous, current;     ///< ray data
        double squared_loss = -1.; ///< squared residual, -1. means it hasn't been evaluated yet

        void set_loss(double r) {
            squared_loss = r * r;
        }

        double get_loss_squared() const {
            return squared_loss;
        }

        // interface to cost function
        std::pair<Ray, Ray> MakePair() const {
            return std::make_pair(previous, current);
        }
    };

    using Matches = std::vector<Match>;
    using MatchesPtr = std::shared_ptr<Matches>;
    using MatchesConstPtr = std::shared_ptr<const Matches>;

    using StorageDataType =
        std::vector<std::tuple<MatchesPtr, Camera::ConstPtr, commons::CostFunctionType>>;

public: // attributes
public: // public methods
    MotionOptimizer();

    /// /////////////////////////////////////////////////////////////////
    /// \brief setter for SolverOptions
    void SetSolverOptions(Options);

    /// /////////////////////////////////////////////////////////////////
    /// \brief getter for SolverOptions
    Options GetSolverOptions();

    /// /////////////////////////////////////////////////////////////////
    /// @brief add data(matches, camera, cost_function_type) to cost functions before solving it
    /// @par matches, non constant ptr to matches, after solving residuals are stored here so you
    /// can access them
    /// @par camera, const_ptr to camera
    /// @par type, type pof the cost functino used in optimization
    /// @todo remove hard coded scaling
    void AddData(MatchesPtr matches, Camera::ConstPtr camera, commons::CostFunctionType type);

    /// /////////////////////////////////////////////
    /// @brief check if optimizer can start to optimize
    /// @return  bool
    /// @par void
    bool IsReadyToOptimize();

    /// /////////////////////////////////////////////////////////////////
    /// @brief delete added data
    void ClearData();

    /// /////////////////////////////////////////////////////////////////
    /// @brief solve the optimization problem
    /// @par x, roll_pitch_yaw whhich will be optimized
    /// @par arc_length, fix_arclength, will not be optimized (@todo make it optimizable in curves)
    /// @par verbose, if true prints stuff to cout
    /// @return returns the ceres::Solver::Summary::FullReport
    std::string Solve(std::vector<double>& x, double arc_length, bool verbose = false);

    /// /////////////////////////////////////////////////////////////////
    /// @brief overload for a client interface who is only intereset in the optimized pose, uses the
    /// last optimized value as start value
    /// @par arc_length, fix_arclength, will not be optimized (@todo make it optimizable in curves)
    /// @par verbose, if true prints stuff to cout
    /// @return optimized motion
    Eigen::Affine3d Solve(double arc_length, bool verbose = false);

    /// /////////////////////////////////////////////////////////////////
    ///@brief define start values for optimization
    //    void SetStartValues(std::vector<double> start_roll_pitch_yaw);

    /// /////////////////////////////////////////////////////////////////
    ///@brief get ceres report(ceres::Solver::Summary::FullReport) of last optimization
    std::string GetReport();

private:
    StorageDataType data_;   ///< storage for ceres data
    Options solver_options_; ///< solver options, same as ceres options but don't want to have
                             /// dependency
    std::vector<double> last_estimate_; ///< start values for optimizaion
    std::string report_;                ///< report of last optimization
    double relative_cost_change_;       ///< store relative cost change to see if solver was able to
                                        /// reduce error if not reset prior
};
}
