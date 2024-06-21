#pragma once
#include "../work_data.hpp"
#include "../configuration.hpp"
#include "../invariants.hpp"

namespace moves {

  class regroup_segment {
    work_data_t &wdata;
    configuration_t &config;
    triqs::mc_tools::random_generator &rng;

    // Internal data
    int color;
    segment_t left_seg, right_seg;
    int left_seg_idx, right_seg_idx;
    bool making_full_line;
    double det_sign;

    public:
    regroup_segment(work_data_t &data_, configuration_t &config_, triqs::mc_tools::random_generator &rng_)
       : wdata(data_), config(config_), rng(rng_) {};
    // ------------------
    double attempt();
    double accept();
    void reject();
  };
}; // namespace moves
