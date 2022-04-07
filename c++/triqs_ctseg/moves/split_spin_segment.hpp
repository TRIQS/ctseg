#pragma once
#include "../work_data.hpp"
#include "../configuration.hpp"

namespace moves {
  class split_spin_segment {
    work_data_t &wdata;
    configuration_t &config;
    triqs::mc_tools::random_generator &rng;

    // Internal dataseg;
    long line_idx, idx_up, idx_down;
    tau_t tau_up, tau_down;
    double ln_trace_ratio, det_ratio, prop_ratio, det_sign;
    std::pair<long, tau_t> propose(int color);

    public:
    split_spin_segment(work_data_t &data_, configuration_t &config_, triqs::mc_tools::random_generator &rng_)
       : wdata(data_), config(config_), rng(rng_) {}

    // ------------------
    double attempt();
    double accept();
    void reject();
  };
}; // namespace moves
