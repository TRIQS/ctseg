#pragma once

#include "./types.hpp"

// Parameters for the solver construction

struct constr_params_t {

  /// Inverse temperature
  double beta;

  /// Structure of the GF (names, sizes of blocks)
  gf_struct_t gf_struct;

  /// Number of time slices for $Delta(\tau)$/$G(\tau)$/$F(\tau)$
  int n_tau = 10001;

  /// Number of time slices for $K(\tau)$
  int n_tau_k = 10001;

  /// Number of time slices for $J_\perp(\tau)$
  int n_tau_jperp = 10001;
};

//---------------------------------------------
// Parameters for the solve method

struct solve_params_t {

  /// local Hamiltonian
  Op h_int;

  /// Number of QMC cycles
  int n_cycles;

  /// Length of a single QMC cycle
  int length_cycle = 50;

  /// Number of cycles for thermalization
  int n_warmup_cycles = 5000;

  /// Seed for random number generator
  int random_seed = 34788 + 928374 * mpi::communicator().rank();

  /// Name of random number generator
  std::string random_name = "";

  /// Maximum runtime in seconds, use -1 to set infinite
  int max_time = -1;

  /// Verbosity level
  int verbosity = mpi::communicator().rank() == 0 ? 3 : 0;

  /// Whether to perform the move insert segment
  bool move_insert_segment_v2 = false;

  /// Whether to perform the move remove segment
  bool move_remove_segment_v2 = false;

  /// Whether to perform the move split segment
  bool move_split_segment_v2 = false;

  /// Whether to perform the move group into spin segment
  bool move_regroup_segment_v2 = false;

  /// Whether to perform the move insert segment
  bool move_insert_segment = true;

  /// Whether to perform the move remove segment
  bool move_remove_segment = true;

  /// Whether to perform the move move segment
  bool move_move_segment = false;

  /// Whether to perform the move split segment
  bool move_split_segment = false;

  /// Whether to perform the move group into spin segment
  bool move_regroup_segment = false;

  /// Whether to perform the move insert spin segment
  bool move_insert_spin_segment = false;

  /// Whether to perform the move remove spin segment
  bool move_remove_spin_segment = false;

  /// Whether to perform the move swap spin lines
  bool move_swap_spin_lines = false;

  /// Whether to measure G(tau) (see [[measure_gt]])
  bool measure_gt = true;

  /// Whether to measure density (see [[measure_density]])
  bool measure_n = true;

  /// Whether to measure <nn> (see [[measure_nn]])
  bool measure_nn = false;

  /// Whether to measure langle n(tau)n(0)rangle (see [[measure_nnt]])
  bool measure_nnt = false;

  /// Hartree shift of the chem pot
  nda::vector<double> hartree_shift = nda::vector<double>{};

  /// The maximum size of the determinant matrix before a resize
  int det_init_size = 100;

  /// Max number of ops before the test of deviation of the det, M^-1 is performed.
  int det_n_operations_before_check = 100;

  /// Threshold for determinant precision warnings
  double det_precision_warning = 1.e-8;

  /// Threshold for determinant precision error
  double det_precision_error = 1.e-5;

  /// Bound for the determinant matrix being singular, abs(det) > singular_threshold. If <0, it is !isnormal(abs(det))
  double det_singular_threshold = -1;

  //solve_params_t(){};
};

/// A struct combining both constr_params_t and solve_params_t
struct params_t : constr_params_t, solve_params_t {
  params_t(constr_params_t const &constr_params_, solve_params_t const &solve_params_)
     : constr_params_t{constr_params_}, solve_params_t{solve_params_} {}
};
