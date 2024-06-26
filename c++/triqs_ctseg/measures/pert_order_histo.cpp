#include "./pert_order_histo.hpp"
#include "../logs.hpp"

namespace triqs_ctseg::measures {

  pert_order_histo::pert_order_histo(params_t const &p, work_data_t const &wdata,
                                                     configuration_t const &config, results_t &results)
     : wdata{wdata}, config{config}, results{results} {

    histo_delta = triqs::stat::histogram{0, p.histogram_max_order};
    histo_Jperp = triqs::stat::histogram{0, p.histogram_max_order};
  }

  // -------------------------------------

  void pert_order_histo::accumulate(double) {

    // For the config, compute the number of segments
    long n_segments =
       std::accumulate(begin(config.seglists), end(config.seglists), 0, [](long r, auto &&v) { return r + v.size(); });
    long n_Splusminus = 2 * config.Jperp_list.size();
    long delta_order =
       n_segments - n_Splusminus; // half # of c, cdag operators  = (2 * n_segments - 2* n_Splusminus) /2

    histo_delta << delta_order;
    histo_Jperp << config.Jperp_list.size();
  }

  // -------------------------------------

  void pert_order_histo::collect_results(mpi::communicator const &c) {

    histo_delta = mpi::all_reduce(histo_delta, c);
    histo_Jperp = mpi::all_reduce(histo_Jperp, c);

    // Normalize
    results.pert_order_histo_Delta = pdf(histo_delta);
    results.pert_order_histo_Jperp = pdf(histo_Jperp);
  }

} // namespace triqs_ctseg::measures
