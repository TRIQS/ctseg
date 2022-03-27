#pragma once
#include "../work_data.hpp"
#include "../configuration.hpp"

namespace moves {

  class move_segment {
    work_data_t &wdata;
    configuration_t &config;
    triqs::mc_tools::random_generator &rng;

    // Internal data
    int origin_color      = 0;
    int destination_color = 0;
    segment_t origin_segment;
    int origin_index{};
    std::vector<segment_t>::const_iterator destination_it;
    qmc_time_factory_t fac = qmc_time_factory_t{wdata.beta};
    double det_sign;
    // Internal methods
    bool do_overlap(segment_t seg1, segment_t seg2);
    bool is_movable(std::vector<segment_t> const &seglist, segment_t const &seg);

    public:
    // Constructor
    move_segment(work_data_t &data_, configuration_t &config_, triqs::mc_tools::random_generator &rng_)
       : wdata(data_), config(config_), rng(rng_){};
    // ------------------
    double attempt();
    double accept();
    void reject();
  };
}; // namespace moves
