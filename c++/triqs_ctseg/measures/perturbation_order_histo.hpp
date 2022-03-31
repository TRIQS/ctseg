#pragma once
#include "../configuration.hpp"
#include "../results.hpp"
#include "../work_data.hpp"

namespace measures {

  struct perturbation_order_histo {

    work_data_t const &wdata;
    configuration_t const &config;
    results_t &results;

    triqs::stat::histogram histo_delta, histo_Jperp;

    perturbation_order_histo(params_t const &params, work_data_t const &wdata, configuration_t const &config, results_t &results);

    void accumulate(double s);
    void collect_results(mpi::communicator const &c);
  };
} // namespace measures
