// Copyright (c) 2022--present, The Simons Foundation
// Copyright (c) 2022--present, Max Planck Institute for Polymer Research, Mainz, Germany
// This file is part of TRIQS/ctseg and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#pragma once

#include <triqs/gfs.hpp>
#include <triqs/operators/many_body_operator.hpp>
using namespace triqs::gfs;

namespace triqs_ctseg {

  /// Parameters used for constructing the solver class.
  struct constr_params_t {

    /// Inverse temperature \f$ \beta \f$.
    double beta;

    /// Structure of the Green's function (names and sizes of blocks).
    gf_struct_t gf_struct;

    /// Number of time slices for fermionic functions.
    int n_tau = 10001;

    /// Number of time slices for bosonic functions.
    int n_tau_bosonic = 10001;
  };

  //---------------------------------------------

  /// Parameters passed to the ``solve()`` method of the solver class.
  struct solve_params_t {

    /// Quartic part of the local Hamiltonian.
    triqs::operators::many_body_operator h_int;

    /// Quandratic part of the local Hamiltonian (including chemical potential).
    triqs::operators::many_body_operator h_loc0;

    /// Number of points on which to measure \f$ G(\tau) \f$ / \f$ F(\tau) \f$ (defaults to ``n_tau``).
    int n_tau_G = 0;

    /// Number of points on which to measure 2-point functions (defaults to ``n_tau_bosonic``.)
    int n_tau_chi2 = 0;

    /// DLR frequency cutoff.
    double dlr_omega_max = 100;

    /// DLR precision.
    double dlr_epsilon = 1e-8;

    /// Number of bosonic M-frequency points on which to measure vertex functions.
    int n_w_b_vertex = 10;

    /// Number of fermionic M-frequency points on which to measure vertex functions.
    int n_w_f_vertex = 10;

    /// Number of QMC cycles.
    int n_cycles;

    /// Length of a single QMC cycle.
    int length_cycle = 50;

    /// Number of cycles for thermalization.
    int n_warmup_cycles = 5000;

    /// Seed for random number generator.
    int random_seed = 34788 + 928374 * mpi::communicator().rank();

    /// Name of random number generator.
    std::string random_name = "";

    /// Maximum runtime in seconds, use -1 to set infinite.
    int max_time = -1;

    /// Verbosity level.
    int verbosity = mpi::communicator().rank() == 0 ? 3 : 0;

    // -------- Move control --------------

    /// Whether to perform the move insert segment.
    bool move_insert_segment = true;

    /// Whether to perform the move remove segment.
    bool move_remove_segment = true;

    /// Whether to perform the move double insert segment.
    bool move_double_insert_segment = true;

    /// Whether to perform the move double remove segment.
    bool move_double_remove_segment = true;

    /// Whether to perform the move move segment.
    bool move_move_segment = true;

    /// Whether to perform the move split segment.
    bool move_split_segment = true;

    /// Whether to perform the move group into spin segment.
    bool move_regroup_segment = true;

    /// Whether to perform the move insert spin segment.
    bool move_insert_spin_segment = true;

    /// Whether to perform the move remove spin segment.
    bool move_remove_spin_segment = true;

    /// Whether to perform the move insert spin segment.
    bool move_split_spin_segment = true;

    /// Whether to perform the move remove spin segment.
    bool move_regroup_spin_segment = true;

    /// Whether to perform the move swap spin lines.
    bool move_swap_spin_lines = true;

    // -------- Measure control --------------

    /// Whether to measure the perturbation order histograms (order in Delta and Jperp).
    bool measure_pert_order = true;

    /// Whether to measure \f$ G(\tau) \f$.
    bool measure_G_tau = true;

    /// Whether to measure \f$ F(\tau) \f$.
    bool measure_F_tau = false;

    /// Whether to measure densities.
    bool measure_densities = true;

    /// Whether to measure the average sign.
    bool measure_average_sign = true;

    /// Whether to measure \f$ \langle n(0) n(0) \rangle \f$.
    bool measure_nn_static = false;

    /// Whether to measure \f$ \langle n(\tau) n(0) \rangle \f$.
    bool measure_nn_tau = false;

    /// Whether to measure \f$ \langle n(\nu)n(0) \rangle \f$.
    bool measure_nn_nu_dlr = false;

    /// Whether to measure \f$ \langle S_x(\tau) S_x(0) \rangle \f$.
    bool measure_Sperp_tau = false;

    /// Whether to measure oriented transverse spin correlations
    /// \f$ \langle S^-(\tau) S^+(0) \rangle \f$ and
    /// \f$ \langle S^+(\tau) S^-(0) \rangle \f$.
    bool measure_Sperp_asym_tau = false;

    /// Whether to measure the occupation-basis TTI diagonal density matrix.
    bool measure_density_matrix = true;

    /// Legacy alias for measure_density_matrix.
    bool measure_state_hist = false;

    /// Whether to measure retarded static correlations for tail moments.
    bool measure_dyn_corr = false;

    /// Whether to measure three-point correlation function.
    bool measure_g2w = false;

    /// Whether to measure four-point correlation function.
    bool measure_g3w = false;

    // -------- Misc parameters --------------

    /// Threshold below which the imaginary part of the local Hamiltonian h_loc0 is set to zero
    /// (CT-SEG uses a real h_loc0); above it the solver errors. Raise to accept a larger
    /// imaginary part.
    double imag_threshold = 1.e-13;

    /// The maximum size of the determinant matrix before a resize.
    int det_init_size = 100;

    /// Max number of ops before testing the accuracy of \f$ \det(M) \f$ and \f$ M^{-1} \f$.
    int det_n_operations_before_check = 100;

    /// Threshold for determinant precision warnings.
    double det_precision_warning = 1.e-8;

    /// Threshold for determinant precision error.
    double det_precision_error = 1.e-5;

    /// Bound for the determinant matrix being singular (if \f$ < 0 \f$, checks for subnormal numbers).
    double det_singular_threshold = -1;

    /// Maximum order for the perturbation order histograms.
    int histogram_max_order = 1000;

    /// Output characteristic configurations in a separate file.
    bool visualize_config = false;
  };

  /// A struct combining both constr_params_t and solve_params_t
  struct params_t : constr_params_t, solve_params_t {
    params_t(constr_params_t const &constr_params_, solve_params_t const &solve_params_)
       : constr_params_t{constr_params_}, solve_params_t{solve_params_} {}
  };

  /// Write all containers to hdf5 file
  void h5_write(h5::group h5group, std::string subgroup_name, constr_params_t const &c);

  /// Reads all containers from hdf5 file
  void h5_read(h5::group h5group, std::string subgroup_name, constr_params_t &c);

  /// Write all containers to hdf5 file
  void h5_write(h5::group h5group, std::string subgroup_name, solve_params_t const &c);

  /// Reads all containers from hdf5 file
  void h5_read(h5::group h5group, std::string subgroup_name, solve_params_t &c);

} // namespace triqs_ctseg
