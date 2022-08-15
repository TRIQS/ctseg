#include "./sperp_tau.hpp"
#include "../logs.hpp"

namespace measures {

  sperp_tau::sperp_tau(params_t const &p, work_data_t const &wdata, configuration_t const &config, results_t &results)
     : wdata{wdata}, config{config}, results{results} {

    beta = p.beta;

    n_color = config.n_color();

    ss_tau   = gf<imtime>({beta, Boson, p.n_tau_k}, {1, 1});
    ss_tau() = 0;
    Z        = 0;
  }

  // -------------------------------------

  void sperp_tau::accumulate(double s) {

    LOG("\n =================== MEASURE < S_x S_x > (tau) ================ \n");

    Z += s;

    for (auto const &[k, line] : itertools::enumerate(config.Jperp_list)) {
      auto dtau1 = double(line.tau_Splus - line.tau_Sminus);
      auto dtau2 = double(line.tau_Sminus - line.tau_Splus);
      ss_tau[closest_mesh_pt(dtau1)] += 0.5 / (real(wdata.Jperp(dtau1)(0, 0)));
      ss_tau[closest_mesh_pt(dtau2)] += 0.5 / (real(wdata.Jperp(dtau2)(0, 0)));
    }

#if 0
    auto N_lines = config.Jperp_list.size();

    for (auto const &[i, line1] : itertools::enumerate(config.Jperp_list)) {
      for (auto const &[j, line2] : itertools::enumerate(config.Jperp_list)) {
        auto Jperp_list_loc         = config.Jperp_list;
        Jperp_list_loc[i].tau_Splus = config.Jperp_list[j].tau_Splus;
        Jperp_list_loc[j].tau_Splus = config.Jperp_list[i].tau_Splus;
        auto t11                    = double(line1.tau_Splus - line1.tau_Sminus);
        auto t22                    = double(line2.tau_Splus - line2.tau_Sminus);
        auto t12                    = double(line1.tau_Splus - line2.tau_Sminus);
        auto t21                    = double(line2.tau_Splus - line1.tau_Sminus);
        auto w_mod                  = real(wdata.Jperp(t12)(0, 0)) * real(wdata.Jperp(t21)(0, 0))
           / (real(wdata.Jperp(t11)(0, 0)) * real(wdata.Jperp(t22)(0, 0)));

        for (auto const &[k, line] : itertools::enumerate(Jperp_list_loc)) {
          auto dtau1 = double(line.tau_Splus - line.tau_Sminus);
          auto dtau2 = double(line.tau_Sminus - line.tau_Splus);
          ss_tau[closest_mesh_pt(dtau1)] += w_mod * 0.5 / (real(wdata.Jperp(dtau1)(0, 0)) * N_lines * N_lines);
          ss_tau[closest_mesh_pt(dtau2)] += w_mod * 0.5 / (real(wdata.Jperp(dtau2)(0, 0)) * N_lines * N_lines);
        }
      }
    }
#endif

  }

  // -------------------------------------

  void sperp_tau::collect_results(mpi::communicator const &c) {

    Z = mpi::all_reduce(Z, c);

    ss_tau = mpi::all_reduce(ss_tau, c);
    ss_tau = ss_tau / (-beta * Z * ss_tau.mesh().delta());

    // Fix the point at zero and beta
    ss_tau[0] *= 2;
    ss_tau[ss_tau.mesh().size() - 1] *= 2;

    // store the result (not reused later, hence we can move it).
    results.sperp_tau = std::move(ss_tau);
  }

} // namespace measures
