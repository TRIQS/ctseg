#pragma once
#include "../configuration.hpp"
#include "../work_data.hpp"
#include "../results.hpp"

namespace measures {

  struct g_f_tau {

    work_data_t const &wdata;
    configuration_t const &config;
    results_t &results;
    double beta;
    bool measure_ft;
    gf_struct_t gf_struct;

    block_gf<imtime> g_tau;
    block_gf<imtime> f_tau;

    double Z;

    g_f_tau(params_t const &params, work_data_t const &wdata, configuration_t const &config, results_t &results);

    void accumulate(double s);
    void collect_results(mpi::communicator const &c);
    double fprefactor(long const &block, std::pair<tau_t, long> const &y);
  };

} // namespace measures
