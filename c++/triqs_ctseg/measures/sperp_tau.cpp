#include "./sperp_tau.hpp"
#include "../logs.hpp"

namespace measures {

  sperp_tau::sperp_tau(params_t const &p, work_data_t const &wdata, configuration_t const &config, results_t &results)
     : wdata{wdata}, config{config}, results{results} {

    beta       = p.beta;

    n_color = config.n_color();

    ss_tau   = gf<imtime>({beta, Boson, p.n_tau_k}, {1,1});
    ss_tau() = 0;
    Z       = 0;
  }

  // -------------------------------------

  void sperp_tau::accumulate(double s) {

    LOG("\n =================== MEASURE < S_x S_x > (tau) ================ \n");

    Z += s;

    for (auto const &[i, line] : itertools::enumerate(config.Jperp_list)) {
        auto dtau = double(line.tau_Splus - line.tau_Sminus);
        ss_tau[closest_mesh_pt(dtau)] += 1.0 / real(wdata.Jperp(dtau)(0,0));

    }  
  }

  // -------------------------------------

  void sperp_tau::collect_results(mpi::communicator const &c) {

    Z = mpi::all_reduce(Z, c);

    ss_tau = mpi::all_reduce(ss_tau, c);
    ss_tau = ss_tau / (- beta * Z * ss_tau.mesh().delta());

    // Fix the point at zero and beta
    ss_tau[0] *= 2;
    ss_tau[ss_tau.mesh().size() - 1] *= 2;

    // store the result (not reused later, hence we can move it).
    results.sperp_tau = std::move(ss_tau);
  }

} // namespace measures
