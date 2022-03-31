#pragma once
#include "../work_data.hpp"
#include "../configuration.hpp"

namespace moves {

  class split_segment {
    work_data_t &wdata;
    configuration_t &config;
    triqs::mc_tools::random_generator &rng;

    // Internal data
    int color;
    segment_t prop_seg;
    dimtime_t tau_left, tau_right;
    int prop_seg_idx, right_seg_idx;
    bool splitting_full_line;
    double det_sign;

    public:
    split_segment(work_data_t &data_, configuration_t &config_, triqs::mc_tools::random_generator &rng_)
       : wdata(data_), config(config_), rng(rng_){};
    // ------------------
    double attempt();
    double accept();
    void reject();
  };
}; // namespace moves
