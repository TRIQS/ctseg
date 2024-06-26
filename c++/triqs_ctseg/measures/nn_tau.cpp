#include "./nn_tau.hpp"
#include "../logs.hpp"

namespace measures {

  nn_tau::nn_tau(params_t const &p, work_data_t const &wdata, configuration_t const &config, results_t &results)
     : wdata{wdata}, config{config}, results{results} {

    beta    = p.beta;
    ntau    = p.n_tau_chi2;
    dtau    = p.beta / (ntau - 1);
    n_color = config.n_color();
    block_number = wdata.block_number;
    index_in_block = wdata.index_in_block;

    q_tau   = gf<imtime>({beta, Boson, ntau}, {n_color, n_color});
    q_tau() = 0;
    q_tau_block = make_block2_gf<imtime>({beta, Boson, ntau}, p.gf_struct);
  }

  // -------------------------------------

  // FIXME mv to nda
  inline static range::all_t nda_all = range::all_t{};

  void nn_tau::accumulate(double s) {

    LOG("\n =================== MEASURE nn(tau) ================ \n");

    Z += s;

    // <n_a(tau) n_b(0) >
    for (int b = 0; b < n_color; ++b) {
      if (n_at_boundary(config.seglists[b]) == 0) continue; // nb = 0, nothing to accumulate
      for (int a = 0; a < n_color; ++a)
        for (auto const &seg : config.seglists[a]) {

          // find closest mesh point to the right of c
          int u_idx_c = (seg.tau_c == tau_t::beta()) ? ntau - 1 : int(std::floor(seg.tau_c / dtau));
          // find closest mesh point to the left of cdag
          int u_idx_cdag = int(std::ceil(seg.tau_cdag / dtau));

          auto d = q_tau.data()(nda_all, a, b); // a view of the data for fixed a,b
          // add + s to the data. NB : id1 > id2
          auto fill = [s, d](long u_idx1, long u_idx2) {
            ALWAYS_EXPECTS((u_idx1 >= u_idx2), "error", 1);
            for (auto u = u_idx1; u >= u_idx2; --u) d(u) += s;
          };

          // Execute with 2 cases : cyclic segment or not
          if (not is_cyclic(seg)) {
            if (u_idx_c >= u_idx_cdag) fill(u_idx_c, u_idx_cdag);
          } else { // cyclic segment
            ALWAYS_EXPECTS((u_idx_cdag >= u_idx_c), "eee", 1);
            fill(u_idx_c, 0);
            fill(ntau - 1, u_idx_cdag);
          }
        }
    }
  }

  // -------------------------------------

  void nn_tau::collect_results(mpi::communicator const &c) {

    Z = mpi::all_reduce(Z, c);

    q_tau = mpi::all_reduce(q_tau, c);
    q_tau = q_tau / Z; //(beta * Z * q_tau.mesh().delta());

    // store the result 
    for (auto const &c1 : range(n_color)) {
      for (auto const &c2 : range(n_color)) {
        q_tau_block(block_number[c1], block_number[c2]).data()(range::all, index_in_block[c1], index_in_block[c2]) = q_tau.data()(range::all, c1, c2);
      }
    }
    results.nn_tau = std::move(q_tau_block);
  }
} // namespace measures
