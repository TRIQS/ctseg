#pragma once
#include "../work_data.hpp"
#include "../configuration.hpp"
#include "../invariants.hpp"

namespace moves {
  class remove_spin_segment {
    work_data_t &wdata;
    configuration_t &config;
    triqs::mc_tools::random_generator &rng;

    // Internal data
    int line_idx, orig_color, dest_color, dest_right_idx, dest_left_idx;
    std::vector<segment_t>::const_iterator orig_it, dest_it;
    segment_t spin_seg;
    bool making_full_line;
    double det_sign;

    public:
    remove_spin_segment(work_data_t &data_, configuration_t &config_, triqs::mc_tools::random_generator &rng_)
       : wdata(data_), config(config_), rng(rng_){};
    // ------------------
    double attempt();
    double accept();
    void reject();
  };
}; // namespace moves
