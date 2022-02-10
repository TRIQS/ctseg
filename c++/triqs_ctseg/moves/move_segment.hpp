#pragma once
#include "../params.hpp"
#include "../qmc_data.hpp"
#include "../configuration.hpp"

namespace moves {

  class move_segment {
    param_t const &params;
    qmc_data const &data;
    configuration &config;
    triqs::mc_tools::random_generator &rng;

    // Internal data // FIXME
    int origin_color = 0;
    int destination_color = 0; 
    segment_t origin_segment;
    std::vector<segment_t>>::iterator origin_index;
    std::vector<segment_t>>::iterator destination_index;

    public:
    // Constructor
    move_segment(const &params_, const qmc_data &data_, configuration &config_, triqs::mc_tools::random_generator &rng_)
      : params(params_), data(data_), config(config_), rng(rng_){};
    // ------------------
    double attempt();
    double accept();
    void reject();
  };
}; // namespace moves