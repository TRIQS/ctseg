#pragma once
#include "../work_data.hpp"
#include "../configuration.hpp"

namespace moves {

  class regroup_segment {
    work_data_t const &wdata;
    configuration_t &config;
    triqs::mc_tools::random_generator &rng;

    // Internal data 
    int color{};
    segment_t left_segment;
    segment_t right_segment;
    int left_segment_index{};
    int right_segment_index{};
    bool making_full_line{};
    time_point_factory_t time_point_factory = time_point_factory_t{wdata.beta};

    public:
    // Constructor
    regroup_segment(const work_data_t &data_, configuration_t &config_, triqs::mc_tools::random_generator &rng_)
      : wdata(data_), config(config_), rng(rng_){};
    // ------------------
    double attempt();
    double accept();
    void reject();
  };
}; // namespace moves