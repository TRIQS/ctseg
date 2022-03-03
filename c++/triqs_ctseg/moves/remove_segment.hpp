#pragma once
#include "../work_data.hpp"
#include "../configuration.hpp"


namespace moves {
  class remove_segment {
    work_data_t const &wdata;
    configuration_t &config;
    triqs::mc_tools::random_generator &rng;

    // Internal data
    int color = 0;
    segment_t proposed_segment;
    int proposed_segment_index{};
    time_point_factory_t time_point_factory = time_point_factory_t{wdata.beta};

    public:
    // Constructor
    remove_segment(const work_data_t &data_, configuration_t &config_, triqs::mc_tools::random_generator &rng_)
      : wdata(data_), config(config_), rng(rng_){};
    // ------------------
    double attempt();
    double accept();
    void reject();
  };
}; // namespace moves