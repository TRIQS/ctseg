#pragma once
#include "../params.hpp"
#include "../qmc_data.hpp"
#include "../configuration.hpp"

using time_point_factory_t = triqs::utility::time_segment;

namespace moves {

  // FIXME : LIGNE PLEINE ?? dans overlap ???

  /// Insert a segment $c_a(\tau)c_a^\dagger(\tau')$
  //
  class insert_segment {
    param_t const &params;
    qmc_data const &data;
    configuration &config;
    triqs::mc_tools::random_generator &rng;

    // internal data

    int color = 0;
    segment_t proposed_segment;
    time_point_factory_t time_point_factory = time_point_factory_t{params.beta};

    public:
    // FIXME
    //insert_segment(const qmc_data *data_, configuration *config_, triqs::mc_tools::random_generator &rng_)
    //  : data(data_), config(config_), rng(rng_){};

    // ----------------------------

    double attempt();
    double accept();
    void reject();
  };
} // namespace moves
