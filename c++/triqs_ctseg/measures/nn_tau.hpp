#pragma once
#include "../configuration.hpp"
#include "../work_data.hpp"
#include "../results.hpp"

namespace triqs_ctseg::measures {

  struct nn_tau {

    work_data_t const &wdata;
    configuration_t const &config;
    results_t &results;
    double beta;
    double dtau;
    int ntau;
    std::vector<long> block_number, index_in_block;

    gf<imtime> q_tau;
    block2_gf<imtime> q_tau_block;

    double Z = 0;
    int n_color;

    nn_tau(params_t const &params, work_data_t const &wdata, configuration_t const &config, results_t &results);

    void accumulate(double s);
    void collect_results(mpi::communicator const &c);
  };

} // namespace triqs_ctseg::measures
