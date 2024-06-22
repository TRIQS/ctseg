#pragma once
#include "../configuration.hpp"
#include "../results.hpp"
#include "../work_data.hpp"

namespace measures {

  struct densities {

    work_data_t const &wdata;
    configuration_t const &config;
    results_t &results;

    nda::array<double, 1> n;

    double Z = 0;

    densities(params_t const &params, work_data_t const &wdata, configuration_t const &config, results_t &results);

    void accumulate(double s);
    void collect_results(mpi::communicator const &c);
  };
} // namespace measures
