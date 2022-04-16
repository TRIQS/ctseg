#include "./g_f_tau.hpp"
#include "../logs.hpp"

namespace measures {

  g_f_tau::g_f_tau(params_t const &p, work_data_t const &wdata, configuration_t const &config, results_t &results)
     : wdata{wdata}, config{config}, results{results} {

    beta       = p.beta;
    measure_ft = p.measure_ft and wdata.rot_inv;
    gf_struct  = p.gf_struct;

    g_tau   = block_gf<imtime>{triqs::mesh::imtime{beta, Fermion, p.n_tau}, p.gf_struct};
    f_tau   = block_gf<imtime>{triqs::mesh::imtime{beta, Fermion, p.n_tau}, p.gf_struct};
    g_tau() = 0;
    f_tau() = 0;
  }

  // -------------------------------------

  void g_f_tau::accumulate(double s) {

    LOG("\n =================== MEASURE G(tau) ================ \n");

    Z += s;

    for (auto [bl_idx, det] : itertools::enumerate(wdata.dets)) {
      foreach (det, [&, &block = bl_idx](auto const &x, auto const &y, auto const &M) {
        auto &g = g_tau[block];
        // beta-periodicity is implicit in the argument, just fix the sign properly
        auto val  = (y.first >= x.first ? s : -s) * M;
        auto dtau = double(y.first - x.first);
        g[closest_mesh_pt(dtau)](y.second, x.second) += val;

        if (measure_ft) {
          auto &f = f_tau[block];
          f[closest_mesh_pt(dtau)](y.second, x.second) += val * fprefactor(block, y);
        }
      })
        ;
    }
  }

  // -------------------------------------

  void g_f_tau::collect_results(mpi::communicator const &c) {

    Z = mpi::all_reduce(Z, c);

    g_tau = mpi::all_reduce(g_tau, c);
    g_tau = g_tau / (-beta * Z * g_tau[0].mesh().delta());

    // Fix the point at zero and beta, for each block
    for (auto &g : g_tau) {
      g[0] *= 2;
      g[g.mesh().size() - 1] *= 2;
    }
    // store the result (not reused later, hence we can move it).
    results.G_tau = std::move(g_tau);

    if (measure_ft) {
      f_tau = mpi::all_reduce(f_tau, c);
      f_tau = f_tau / (-beta * Z * f_tau[0].mesh().delta());

      for (auto &f : f_tau) {
        f[0] *= 2;
        f[f.mesh().size() - 1] *= 2;
      }
      results.F_tau = std::move(f_tau);
    }
  }

  double g_f_tau::fprefactor(long const &block, std::pair<tau_t, long> const &y) {
    int color    = wdata.block_to_color(block, y.second);
    double I_tau = 0;
    for (auto const &[c, sl] : itertools::enumerate(config.seglists)) {
      auto ntau = n_tau(y.first, sl); // Density to the right of y.first in sl
      if (c != color) I_tau += wdata.U(c, color) * ntau;
      if (wdata.has_Dt) {
        I_tau -= K_overlap(sl, y.first, false, wdata.Kprime, c, color);
        if (c == color) I_tau -= 2 * real(wdata.Kprime(0)(c, c));
      }
      if (wdata.has_jperp) {
        I_tau -= 4 * real(wdata.Kprime_spin(0)(c, color)) * ntau;
        I_tau -= 2 * K_overlap(sl, y.first, false, wdata.Kprime_spin, c, color);
      }
    }
    return I_tau;
  }
} // namespace measures
