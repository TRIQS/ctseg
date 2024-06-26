#pragma once
#include "../work_data.hpp"
#include "../configuration.hpp"
#include "../invariants.hpp"

namespace triqs_ctseg::moves {

  class insert_spin_segment {
    work_data_t &wdata;
    configuration_t &config;
    triqs::mc_tools::random_generator &rng;

    // Internal data
    int orig_color, dest_color;
    segment_t prop_seg, spin_seg;
    int prop_seg_idx;
    bool splitting_full_line;
    tau_t tau_left, tau_right;
    double det_sign;

    public:
    insert_spin_segment(work_data_t &data_, configuration_t &config_, triqs::mc_tools::random_generator &rng_);
    double attempt();
    double accept();
    void reject();
  };

} // namespace triqs_ctseg::moves
