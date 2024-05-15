#pragma once
#include "../configuration.hpp"
#include "../results.hpp"
#include "../work_data.hpp"

namespace measures {

  struct density {

    work_data_t const &wdata;
    configuration_t const &config;
    results_t &results;

    nda::array<double, 1> densities;

    double Z = 0;

    density(params_t const &params, work_data_t const &wdata, configuration_t const &config, results_t &results);

    void accumulate(double s);
    void collect_results(mpi::communicator const &c);
  };
} // namespace measures
