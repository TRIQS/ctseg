#include "density.hpp"
#include <itertools/itertools.hpp>
#include "../logs.hpp"

namespace measures {

  density::density(params_t const &, work_data_t const &wdata, configuration_t const &config, results_t &results)
     : wdata{wdata}, config{config}, results{results} {

    densities = nda::zeros<double>(config.n_color());
  }

  // -------------------------------------

  void density::accumulate(double s) {

    Z += s;
    for (auto const &[c, seglist] : itertools::enumerate(config.seglists)) {
      double sum = 0;
      for (auto &seg : seglist) sum += double(seg.length()); // accounts for cyclicity
      densities[c] += s * sum;
    }
  }

  // -------------------------------------

  void density::collect_results(mpi::communicator const &c) {
    Z         = mpi::all_reduce(Z, c);
    densities = mpi::all_reduce(densities, c);
    densities /= (Z * tau_t::beta());
    results.densities = densities;
    if (c.rank() == 0) SPDLOG_INFO("Density {}", densities);
  }
} // namespace measures
