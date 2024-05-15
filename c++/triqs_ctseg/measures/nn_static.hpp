#pragma once
#include "../configuration.hpp"
#include "../work_data.hpp"
#include "../results.hpp"

namespace measures {

  struct nn_static {

    work_data_t const &wdata;
    configuration_t const &config;
    results_t &results;
    double beta;

    nda::matrix<double> nn;

    double Z = 0;
    int n_color;

    nn_static(params_t const &params, work_data_t const &wdata, configuration_t const &config, results_t &results);

    void accumulate(double s);
    void collect_results(mpi::communicator const &c);
  };
} // namespace measures
