#pragma once
#include "../params.hpp"
#include "../qmc_data.hpp"
#include "../configuration.hpp"


namespace moves {
  // FIXME : LIGNE PLEINE ?? dans overlap ???
  class insert_segment {
    param_t const &params;
    qmc_data const &data;
    configuration &config;
    triqs::mc_tools::random_generator &rng;

    // Internal data
    int color = 0;
    segment_t proposed_segment;
    std::vector<segment_t>>::iterator proposed_segment_insert_pos;
    time_point_factory_t time_point_factory = time_point_factory_t{params.beta};

    public:
    // Constructor
    insert_segment(const qmc_data *data_, configuration *config_, triqs::mc_tools::random_generator &rng_)
      : data(data_), config(config_), rng(rng_){};
    // ------------------
    double attempt();
    double accept();
    void reject();
  };
}; // namespace moves
