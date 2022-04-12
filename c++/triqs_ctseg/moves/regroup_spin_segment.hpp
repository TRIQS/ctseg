#pragma once
#include "../work_data.hpp"
#include "../configuration.hpp"
#include "../invariants.hpp"

namespace moves {
  class regroup_spin_segment {
    work_data_t &wdata;
    configuration_t &config;
    triqs::mc_tools::random_generator &rng;

    // Internal data
    long idx_c_up, idx_cdag_down, idx_c_down, idx_cdag_up;
    tau_t tau_up, tau_down;
    double ln_trace_ratio, prop_ratio, det_sign;
    bool prop_failed;
    std::tuple<long, long, tau_t, bool> propose(int color);

    public:
    regroup_spin_segment(work_data_t &data_, configuration_t &config_, triqs::mc_tools::random_generator &rng_)
       : wdata(data_), config(config_), rng(rng_) {}

    // ------------------
    double attempt();
    double accept();
    void reject();
  };
}; // namespace moves
