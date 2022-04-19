#include "./nn_tau.hpp"
#include "../logs.hpp"

// FIXME : move this to TRIQS ...
namespace triqs::mesh {

  template <typename Mesh, typename T> typename Mesh::index_t closest_mesh_pt_index(Mesh const &m, T const &x) {
    return closest_point<Mesh, scalar_valued>::invoke(m, closest_mesh_pt(x));
  }

} // namespace triqs::mesh

namespace measures {

  nn_tau::nn_tau(params_t const &p, work_data_t const &wdata, configuration_t const &config, results_t &results)
     : wdata{wdata}, config{config}, results{results} {

    beta    = p.beta;
    n_color = config.n_color();

    q_tau   = gf<imtime>({beta, Boson, p.n_tau}, {n_color, n_color});
    q_tau() = 0;
  }

  // -------------------------------------

  // FIXME mv to nda
  inline static range::all_t nda_all = range::all_t{};

  void nn_tau::accumulate(double s) {

    LOG("\n =================== MEASURE NN(tau) ================ \n");

    Z += s;

    // <n_a(tau) n_b(0) >
    for (int b = 0; b < n_color; ++b) {
      if (n_at_boundary(config, b) == 0) continue; // nb = 0, nothing to accumulate
      for (int a = 0; a < n_color; ++a)
        for (auto const &seg : config.seglists[a]) {

          // position of c cdag of the segment, as an integer index in the mesh
          auto u_idx_c    = closest_mesh_pt_index(q_tau.mesh(), seg.tau_c);
          auto u_idx_cdag = closest_mesh_pt_index(q_tau.mesh(), seg.tau_cdag);

          auto d = q_tau.data()(nda_all, a, b); // a view of the data for fixed a,b
          // add + s to the data. NB : id1 > id2
          auto fill = [s, d](long u_idx1, long u_idx2) {
            ALWAYS_EXPECTS((u_idx1 >= u_idx2), "error", 1);
            for (auto u = u_idx1; u >= u_idx2; --u) d(u) += s;
          };

          // Execute with 2 cases : cyclic segment or not
          if (not is_cyclic(seg)) {
            fill(u_idx_c, u_idx_cdag);
          } else { // cyclic segment
            ALWAYS_EXPECTS((u_idx_cdag >= u_idx_c), "eee", 1);
            fill(u_idx_c, 0);
            fill(q_tau.mesh().size() - 1, u_idx_cdag);
          }
        }
    }
  }

  // -------------------------------------

  void nn_tau::collect_results(mpi::communicator const &c) {

    Z = mpi::all_reduce(Z, c);

    q_tau = mpi::all_reduce(q_tau, c);
    q_tau = q_tau / Z; //(beta * Z * q_tau.mesh().delta());

    // FIXME : do I need this ??
    // Fix the point at zero and beta, for each block
    //q_tau[0] *= 2;
    //q_tau[q_tau.mesh().size() - 1] *= 2;

    // store the result (not reused later, hence we can move it).
    results.nn_tau = std::move(q_tau);
  }
} // namespace measures
