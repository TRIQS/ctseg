#include "./nn_tau.hpp"

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

    q_tau   = gf<imtime>({beta, Boson, p.n_tau}, {wdata.n_color, wdata.n_color});
    q_tau() = 0;
  }

  // -------------------------------------

  void nn_tau::accumulate(double s) {
    Z += s;

    // I take the whole configuration and make a time ordered list of operators.
    auto t_ordered_op_list = make_time_ordered_op_list(config);
    auto zero = time_point_factory.get_lower_pt();
    t_ordered_op_list.emplace_back(zero, 0, false); // add the point 0 to make it easier later

    // the state at 0 or beta
    auto state_at_0   = boundary_state(config);
    auto state_at_tau = state_at_0;
    nda::matrix<double> nn_mat(n_color, n_color);

    long idx1 = q_tau.mesh().size();

    for (auto const &[tau2, color, dagger] : t_ordered_op_list) {
      long idx2 = closest_mesh_pt_index(this->q_tau.mesh(), tau2);

      for (int a = 0; a < n_color; ++a)
        for (int b = 0; b < n_color; ++b) nn_mat(a, b) += int(state_at_tau[a]) * int(state_at_0[b]);

      for (long u = idx1 - 1; u >= idx2; --u) q_tau.data()(u, nda::ellipsis{}) += nn_mat;

      state_at_tau[color] = not state_at_tau[color]; // cross a C or a C dagger ...
      idx1                = idx2;
      ALWAYS_EXPECTS(idx2 < idx1, "Error at {}", tau2);
    }
  }

  // -------------------------------------

  void nn_tau::collect_results(mpi::communicator const &c) {

    Z = mpi::all_reduce(Z, c);

    q_tau = mpi::all_reduce(q_tau, c);
    q_tau = q_tau / (-beta * Z * q_tau.mesh().delta());

    // FIXME : do I need this ??
    // Fix the point at zero and beta, for each block
    //q_tau[0] *= 2;
    //q_tau[q_tau.mesh().size() - 1] *= 2;

    // store the result (not reused later, hence we can move it).
    results.nn_tau = std::move(q_tau);
  }
} // namespace measures
