#pragma once
#include "../params.hpp"
#include "../qmc_data.hpp"
#include "../configuration.hpp"

namespace moves {

  class regroup_segment {
    param_t const &params;
    qmc_data const &data;
    configuration &config;
    triqs::mc_tools::random_generator &rng;

    // Internal data 
    int color = 0;
    segment_t left_segment;
    segment_t right_segment;
    std::vector<segment_t>>::iterator left_segment_index;
    std::vector<segment_t>>::iterator right_segment_index;
    time_point_factory_t time_point_factory = time_point_factory_t{params.beta};

    public:
    // Constructor
    regroup_segment(const &params_, const qmc_data &data_, configuration &config_, triqs::mc_tools::random_generator &rng_)
      : params(params_), data(data_), config(config_), rng(rng_){};
    // ------------------
    double attempt();
    double accept();
    void reject();
  };
}; // namespace moves