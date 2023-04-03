#pragma once
#include "../configuration.hpp"
#include "../results.hpp"
#include "../work_data.hpp"

namespace measures {

  struct sign {

    work_data_t const &wdata;
    configuration_t const &config;
    results_t &results;

    double N = 0; 
    double Z = 0;

    sign(params_t const &params, work_data_t const &wdata, configuration_t const &config, results_t &results);

    void accumulate(double s);
    void collect_results(mpi::communicator const &c);
  };
} // namespace measures
