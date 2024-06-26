#pragma once
#include "../configuration.hpp"
#include "../work_data.hpp"
#include "../results.hpp"
//#include "../precompute_fprefactor.hpp"

namespace triqs_ctseg::measures {

  struct sperp_tau {

    work_data_t const &wdata;
    configuration_t const &config;
    results_t &results;
    double beta;

    gf<imtime> ss_tau;

    double Z = 0;
    int n_color;

    sperp_tau(params_t const &params, work_data_t const &wdata, configuration_t const &config, results_t &results);

    void accumulate(double s);
    void collect_results(mpi::communicator const &c);
  };

} // namespace triqs_ctseg::measures
