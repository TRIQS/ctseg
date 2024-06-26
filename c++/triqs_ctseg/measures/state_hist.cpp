#include "./state_hist.hpp"
#include "../logs.hpp"

namespace triqs_ctseg::measures {

  state_hist::state_hist(params_t const &p, work_data_t const &wdata, configuration_t const &config, results_t &results)
     : wdata{wdata}, config{config}, results{results} {

    beta = p.beta;
    H    = nda::zeros<double>(ipow(2, config.n_color()));
  }

  // -------------------------------------

  void state_hist::accumulate(double s) {

    LOG("\n =================== MEASURE STATE HISTOGRAM  ================ \n");

    /// Measure the state histogram
    /**
    *Measures the time the impurity spends in a certain atomic state:
    * an atomic state is characterized by occupation number for each color, e.g.
    * $|n_1 n_2 n_3 n_4\rangle = |0 1 1 0\rangle$ for a model with 4 colors
    *
    * - the index of a state in the histogram is given by $\sum_i n_i 2^i$
    *
    * - the length of the histogram is 2^n_colors
    */

    Z += s;

    double tau_prev         = beta; // time of prevous operator; start with beta
    nda::vector<bool> state = nda::zeros<bool>(config.n_color());
    for (auto const &op : colored_ordered_ops(config.seglists)) {
      int state_idx = 0;
      for (auto c : range(config.n_color()))
        if (state(c)) state_idx += ipow(2, c); // get the index of the impurity state
      H(state_idx) += (tau_prev - op.tau);
      tau_prev = (double)op.tau;
      ALWAYS_EXPECTS((state(op.color) == op.is_cdag), "Operator error at color {}", op.color);
      state(op.color) = !op.is_cdag;
    }

    // get edge state contribution; tau_prev has time of last operator
    ALWAYS_EXPECTS((state == nda::zeros<bool>(config.n_color())), "Operator error");
    H(0) += tau_prev;
  }
  // -------------------------------------

  void state_hist::collect_results(mpi::communicator const &c) {

    Z = mpi::all_reduce(Z, c);
    H = mpi::all_reduce(H, c);
    H = H / (Z * beta);

    // store the result (not reused later, hence we can move it).
    results.state_hist = std::move(H);
  }

} // namespace triqs_ctseg::measures
