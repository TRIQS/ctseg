#pragma once
#include "../work_data.hpp"
#include "../configuration.hpp"
#include "../invariants.hpp"

namespace moves {
  class remove_segment {
    work_data_t &wdata;
    configuration_t &config;
    triqs::mc_tools::random_generator &rng;

    // Internal data
    int color = 0;
    segment_t prop_seg;
    int prop_seg_idx;
    double det_sign;

    public:
    remove_segment(work_data_t &data_, configuration_t &config_, triqs::mc_tools::random_generator &rng_)
       : wdata(data_), config(config_), rng(rng_){};
    // ------------------
    double attempt();
    double accept();
    void reject();
  };
}; // namespace moves
