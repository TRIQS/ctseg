#include "sign.hpp"
#include <itertools/itertools.hpp>
#include "../logs.hpp"

namespace measures {

  sign::sign(params_t const &, work_data_t const &wdata, configuration_t const &config, results_t &results)
     : wdata{wdata}, config{config}, results{results} {
    Z = 0.0;
    N = 0.0;
  }

  // -------------------------------------

  void sign::accumulate(double s) {
    Z += s;
    N += 1.0;
  }

  // -------------------------------------

  void sign::collect_results(mpi::communicator const &c) {
    Z            = mpi::all_reduce(Z, c);
    N            = mpi::all_reduce(N, c);
    results.sign = Z / N;
  }
} // namespace measures