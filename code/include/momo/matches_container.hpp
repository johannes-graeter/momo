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
//standard includes
#include <string>
#include <vector>
#include <map>
#include <tuple>
#include <memory>

/**
*  @class MatchesContainer
*
*  Stores matches and does thresholding and other operations on them
*
*  Uses google c++ style guide naming conventions
*/
class MatchesContainer {
public: // public classes/enums/types etc...
    using TrackletId = size_t;
    using CameraPoint = std::pair<double, double>;
    using Tracklet = std::vector<CameraPoint>;
    using Matches = std::vector<Tracklet>;
    using MatchesMap = std::map<TrackletId, Tracklet>;

private: // attributes
    MatchesMap matches_map_;
    std::vector<TrackletId> valid_match_indices_;
    bool unset_;

public: // public methods
    // default constructor
    MatchesContainer() : unset_(true) {
        ;
    }

    ///@brief basic setter for matches, set valid indices so that without thresholding,
    /// GetMatches==GetValidMatches
    void SetMatches(const MatchesMap& m) {
        matches_map_ = m;
        unset_ = false;

        SetIndicesAllValid();
    }

    ///@brief store matches in a mapset valid indices so that without thresholding,
    /// GetMatches==GetValidMatches
    void SetMatches(const Matches& m) {
        // if there are no unique ids given from outside assign proper ones
        for (size_t j = 0; j < m.size(); ++j) {
            matches_map_[j] = m[j];
        }

        unset_ = false;

        SetIndicesAllValid();
    }

    ///@brief thresholds the matches by costs and threshold, costs must be in same order as input
    /// matches
    void ThresholdMatches(const Eigen::VectorXd& costs, double threshold) {
        Assert(matches_map_.size() == 0 || int(costs.rows()) == int(matches_map_.size()),
               "thresholding conditions not met");
        valid_match_indices_.clear();
        // preprocess();
        size_t ind = 0;
        for (const auto& el : matches_map_) {
            if (costs[ind] < threshold)
                valid_match_indices_.push_back(el.first);
            ind++;
        }
        valid_match_indices_.shrink_to_fit();
    }

    ///@brief get thresholded matches
    Matches GetMatches() const {
        Assert(!unset_, "in GetMatches - unset");
        Matches cur_m;

        std::transform(matches_map_.cbegin(), matches_map_.cend(), std::back_inserter(cur_m),
                       [](const auto& el) { return el.second; });

        return std::move(cur_m);
    }

    ///@brief get thresholded matches
    Matches GetValidMatches() const {
        Assert(!unset_, "in GetValidMatches - unset");
        Matches cur_m;
        cur_m.reserve(valid_match_indices_.size());
        for (const auto& ind : valid_match_indices_) {
            cur_m.push_back(matches_map_.at(ind));
        }
        // std::transform(valid_match_indices.cbegin(), valid_match_indices.cend(),
        // std::back_inserter(cur_m), [&matches_map_] (const auto& ind){return matches_map_[ind];});
        return std::move(cur_m);
    }

    ///@brief get ids of input matches that are bekow threshold
    std::vector<size_t> GetValidMatchIds() const {
        return valid_match_indices_;
    }

private:
    ///@brief set all matche indices valid
    void SetIndicesAllValid() {
        valid_match_indices_.clear();
        for (size_t j = 0; j < m.size(); ++j) {
            valid_match_indices_.push_back(j);
        }
        valid_match_indices_.shrink_to_fit();
    }

    ///@brief defines custom assert with exception. You can overload it with assert or custom
    ///exception if you want.
    void Assert(bool condition, std::string message) {
        if (!condition)
            throw std::runtime_error("Assertion in MatchesContainer: " + message);
    }
}; // end of MatchesContainer
