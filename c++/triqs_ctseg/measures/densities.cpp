#include "densities.hpp"
#include <itertools/itertools.hpp>
#include "../logs.hpp"

namespace triqs_ctseg::measures {

  densities::densities(params_t const &, work_data_t const &wdata, configuration_t const &config, results_t &results)
     : wdata{wdata}, config{config}, results{results} {

    n = nda::zeros<double>(config.n_color());
  }

  // -------------------------------------

  void densities::accumulate(double s) {

    Z += s;
    for (auto const &[c, seglist] : itertools::enumerate(config.seglists)) {
      double sum = 0;
      for (auto &seg : seglist) sum += double(seg.length()); // accounts for cyclicity
      n[c] += s * sum;
    }
  }

  // -------------------------------------

  void densities::collect_results(mpi::communicator const &c) {
    Z = mpi::all_reduce(Z, c);
    n = mpi::all_reduce(n, c);
    n /= (Z * tau_t::beta());
    results.densities = n;
    if (c.rank() == 0) SPDLOG_INFO("Densities: {}", n);
  }

} // namespace triqs_ctseg::measures
