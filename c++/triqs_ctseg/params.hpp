// Copyright (c) 2022-2024 Simons Foundation
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You may obtain a copy of the License at
//     https://www.gnu.org/licenses/gpl-3.0.txt
//
// Authors: Nikita Kavokine, Hao Lu, Olivier Parcollet, Nils Wentzell

#pragma once

#include <triqs/gfs.hpp>
#include <triqs/operators/many_body_operator.hpp>
using namespace triqs::gfs;

namespace triqs_ctseg {

  // Parameters for the solver construction

  struct constr_params_t {

    /// Inverse temperature
    double beta;

    /// Structure of the Green's function (names and sizes of blocks)
    gf_struct_t gf_struct;

    /// Number of time slices for fermionic functions
    int n_tau = 10001;

    /// Number of time slices for bosonic functions
    int n_tau_bosonic = 10001;
  };

  //---------------------------------------------
  // Parameters for the solve method

  struct solve_params_t {

    /// Quartic part of the local Hamiltonian
    triqs::operators::many_body_operator h_int;

    /// Quandratic part of the local Hamiltonian (including chemical potential)
    triqs::operators::many_body_operator h_loc0;

    /// Number of points on which to measure G(tau)/F(tau) (defaults to n_tau)
    int n_tau_G = 0;

    /// Number of points on which to measure 2-point functions (defaults to n_tau_bosonic)
    int n_tau_chi2 = 0;

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

    /// Whether to measure the perturbation order histograms (order in Delta and Jperp)
    bool measure_pert_order = true;

    /// Whether to measure G(tau) (see measures/g_f_tau)
    bool measure_G_tau = true;

    /// Whether to measure F(tau) (see measures/g_f_tau)
    bool measure_F_tau = false;

    /// Whether to measure densities (see measures/densities)
    bool measure_densities = true;

    /// Whether to measure sign (see measures/sign)
    bool measure_average_sign = true;

    /// Whether to measure <n(0)n(0)> (see measures/nn_static)
    bool measure_nn_static = false;

    /// Whether to measure <n(tau)n(0)> (see measures/nn_tau)
    bool measure_nn_tau = false;

    /// Whether to measure <s_x(tau)s_x(0)> (see measures/sperp_tau)
    bool measure_sperp_tau = false;

    /// Whether to measure state histograms (see measures/state_hist)
    bool measure_state_hist = false;

    // -------- Misc parameters --------------

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

    /// Maximum order for the perturbation order histograms
    int histogram_max_order = 1000;
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
