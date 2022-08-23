#pragma once
#include "../configuration.hpp"
#include "../work_data.hpp"
#include "../results.hpp"
//#include "../precompute_fprefactor.hpp"

namespace measures {

  struct sperp_tau2 {

    work_data_t const &wdata;
    configuration_t const &config;
    results_t &results;
    double beta;

    gf<imtime> ss_tau2;
    gf<imtime> ss_tau2_rev;

    double Z = 0;
    int n_color;

    sperp_tau2(params_t const &params, work_data_t const &wdata, configuration_t const &config, results_t &results);

    void accumulate(double s);
    void collect_results(mpi::communicator const &c);
  };
} // namespace measures