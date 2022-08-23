#pragma once

#include <triqs/gfs.hpp>
#include <triqs/operators/many_body_operator.hpp>
using namespace triqs::gfs;

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
};

//---------------------------------------------
// Parameters for the solve method

struct solve_params_t {

  /// local Hamiltonian
  triqs::operators::many_body_operator h_int;

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

  // -------- Move control --------------

  /// Whether to perform the move insert segment
  bool move_insert_segment = true;

  /// Whether to perform the move remove segment
  bool move_remove_segment = true;

  /// Whether to perform the move move segment
  bool move_move_segment = true;

  /// Whether to perform the move split segment
  bool move_split_segment = true;

  /// Whether to perform the move group into spin segment
  bool move_regroup_segment = true;

  /// Whether to perform the move insert spin segment
  bool move_insert_spin_segment = true;

  /// Whether to perform the move remove spin segment
  bool move_remove_spin_segment = true;

  /// Whether to perform the move insert spin segment
  bool move_split_spin_segment = true;

  /// Whether to perform the move remove spin segment
  bool move_regroup_spin_segment = true;

  /// Whether to perform the move swap spin lines
  bool move_swap_spin_lines = true;

  // -------- Measure control --------------

  /// Whether to measure the perturbation order histograms (Order in Delta, and Jperp)
  bool measure_perturbation_order_histograms = true;

  /// Whether to measure G(tau) (see [[measure_g_f_tau]])
  bool measure_gt = true;

  /// Whether to measure F(tau) (see [[measure_g_f_tau]])
  bool measure_ft = false;

  /// Whether to measure density (see [[measure_density]])
  bool measure_n = true;

  /// Whether to measure <nn> (see [[measure_nn]])
  bool measure_nn = false;

  /// Whether to measure langle n(tau)n(0)rangle (see [[measure_nnt]])
  bool measure_nnt = false;

  /// Whether to measure langle s_x(tau)s_x(0)rangle (see [[measure_sperp_tau]])
  bool measure_sperpt = false;

  /// Whether to measure langle s_x(tau)s_x(0)rangle using N^2 measurements (see [[measure_sperp_tau2]])
  bool measure_sperpt2 = false;

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
};

/// A struct combining both constr_params_t and solve_params_t
struct params_t : constr_params_t, solve_params_t {
  params_t(constr_params_t const &constr_params_, solve_params_t const &solve_params_)
     : constr_params_t{constr_params_}, solve_params_t{solve_params_} {}
};

/// Get the number of colors from the gf_struct.
inline int count_colors(gf_struct_t const &gf_struct) {
  int n = 0;
  for (auto const &[bl_name, bl_size] : gf_struct) { n += bl_size; }
  return n;
};

/// writes all containers to hdf5 file
void h5_write(h5::group h5group, std::string subgroup_name, constr_params_t const &c);

/// reads all containers to hdf5 file
void h5_read(h5::group h5group, std::string subgroup_name, constr_params_t &c);
