// Copyright (c) 2022--present, The Simons Foundation
// Copyright (c) 2022--present, Max Planck Institute for Polymer Research, Mainz, Germany
// This file is part of TRIQS/ctseg and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#include "params.hpp"

namespace triqs_ctseg {

  void h5_write(h5::group h5group, std::string subgroup_name, constr_params_t const &c) {

    h5::group grp = subgroup_name.empty() ? h5group : h5group.create_group(subgroup_name);

    h5_write(grp, "beta", c.beta);
    h5_write(grp, "gf_struct", c.gf_struct);
    h5_write(grp, "n_tau", c.n_tau);
    h5_write(grp, "n_tau_bosonic", c.n_tau_bosonic);
    h5_write(grp, "n_l", c.n_l);
  }

  //------------------------------------

  void h5_read(h5::group h5group, std::string subgroup_name, constr_params_t &c) {

    h5::group grp = subgroup_name.empty() ? h5group : h5group.open_group(subgroup_name);

    h5_read(grp, "beta", c.beta);
    h5_read(grp, "gf_struct", c.gf_struct);
    h5_read(grp, "n_tau", c.n_tau);
    h5_read(grp, "n_tau_bosonic", c.n_tau_bosonic);
    h5::try_read(grp, "n_l", c.n_l);
  }

  //------------------------------------

  void h5_write(h5::group h5group, std::string subgroup_name, solve_params_t const &c) {

    h5::group grp = subgroup_name.empty() ? h5group : h5group.create_group(subgroup_name);

    h5_write(grp, "h_int", c.h_int);
    h5_write(grp, "h_loc0", c.h_loc0);
    h5_write(grp, "n_tau_G", c.n_tau_G);
    h5_write(grp, "n_tau_chi2", c.n_tau_chi2);
    h5_write(grp, "dlr_omega_max", c.dlr_omega_max);
    h5_write(grp, "dlr_epsilon", c.dlr_epsilon);
    h5_write(grp, "n_w_b_vertex", c.n_w_b_vertex);
    h5_write(grp, "n_w_f_vertex", c.n_w_f_vertex);
    h5_write(grp, "n_cycles", c.n_cycles);
    h5_write(grp, "length_cycle", c.length_cycle);
    h5_write(grp, "n_warmup_cycles", c.n_warmup_cycles);
    h5_write(grp, "random_seed", c.random_seed);
    h5_write(grp, "random_name", c.random_name);
    h5_write(grp, "max_time", c.max_time);
    h5_write(grp, "verbosity", c.verbosity);
    h5_write(grp, "move_insert_segment", c.move_insert_segment);
    h5_write(grp, "move_remove_segment", c.move_remove_segment);
    h5_write(grp, "move_double_insert_segment", c.move_double_insert_segment);
    h5_write(grp, "move_double_remove_segment", c.move_double_remove_segment);
    h5_write(grp, "move_move_segment", c.move_move_segment);
    h5_write(grp, "move_split_segment", c.move_split_segment);
    h5_write(grp, "move_regroup_segment", c.move_regroup_segment);
    h5_write(grp, "move_insert_spin_segment", c.move_insert_spin_segment);
    h5_write(grp, "move_remove_spin_segment", c.move_remove_spin_segment);
    h5_write(grp, "move_split_spin_segment", c.move_split_spin_segment);
    h5_write(grp, "move_regroup_spin_segment", c.move_regroup_spin_segment);
    h5_write(grp, "move_swap_spin_lines", c.move_swap_spin_lines);
    h5_write(grp, "measure_pert_order", c.measure_pert_order);
    h5_write(grp, "measure_G_tau", c.measure_G_tau);
    h5_write(grp, "measure_F_tau", c.measure_F_tau);
    h5_write(grp, "measure_G_l", c.measure_G_l);
    h5_write(grp, "measure_F_l", c.measure_F_l);
    h5_write(grp, "measure_densities", c.measure_densities);
    h5_write(grp, "measure_average_sign", c.measure_average_sign);
    h5_write(grp, "measure_nn_static", c.measure_nn_static);
    h5_write(grp, "measure_nn_tau", c.measure_nn_tau);
    h5_write(grp, "measure_nn_nu_dlr", c.measure_nn_nu_dlr);
    h5_write(grp, "measure_Sperp_tau", c.measure_Sperp_tau);
    h5_write(grp, "measure_state_hist", c.measure_state_hist);
    h5_write(grp, "measure_g2w", c.measure_g2w);
    h5_write(grp, "measure_g3w", c.measure_g3w);
    h5_write(grp, "imag_threshold", c.imag_threshold);
    h5_write(grp, "det_init_size", c.det_init_size);
    h5_write(grp, "det_n_operations_before_check", c.det_n_operations_before_check);
    h5_write(grp, "det_precision_warning", c.det_precision_warning);
    h5_write(grp, "det_precision_error", c.det_precision_error);
    h5_write(grp, "det_singular_threshold", c.det_singular_threshold);
    h5_write(grp, "histogram_max_order", c.histogram_max_order);
    h5_write(grp, "visualize_config", c.visualize_config);
  }

  //------------------------------------

  void h5_read(h5::group h5group, std::string subgroup_name, solve_params_t &c) {

    h5::group grp = subgroup_name.empty() ? h5group : h5group.open_group(subgroup_name);

    h5_read(grp, "h_int", c.h_int);
    h5_read(grp, "h_loc0", c.h_loc0);
    h5_read(grp, "n_tau_G", c.n_tau_G);
    h5_read(grp, "n_tau_chi2", c.n_tau_chi2);
    h5_read(grp, "dlr_omega_max", c.dlr_omega_max);
    h5_read(grp, "dlr_epsilon", c.dlr_epsilon);
    h5_read(grp, "n_cycles", c.n_cycles);
    h5_read(grp, "length_cycle", c.length_cycle);
    h5_read(grp, "n_warmup_cycles", c.n_warmup_cycles);
    h5_read(grp, "random_seed", c.random_seed);
    h5_read(grp, "random_name", c.random_name);
    h5_read(grp, "max_time", c.max_time);
    h5_read(grp, "verbosity", c.verbosity);
    h5_read(grp, "move_insert_segment", c.move_insert_segment);
    h5_read(grp, "move_remove_segment", c.move_remove_segment);
    h5_read(grp, "move_double_insert_segment", c.move_double_insert_segment);
    h5_read(grp, "move_double_remove_segment", c.move_double_remove_segment);
    h5_read(grp, "move_move_segment", c.move_move_segment);
    h5_read(grp, "move_split_segment", c.move_split_segment);
    h5_read(grp, "move_regroup_segment", c.move_regroup_segment);
    h5_read(grp, "move_insert_spin_segment", c.move_insert_spin_segment);
    h5_read(grp, "move_remove_spin_segment", c.move_remove_spin_segment);
    h5_read(grp, "move_split_spin_segment", c.move_split_spin_segment);
    h5_read(grp, "move_regroup_spin_segment", c.move_regroup_spin_segment);
    h5_read(grp, "move_swap_spin_lines", c.move_swap_spin_lines);
    h5_read(grp, "measure_pert_order", c.measure_pert_order);
    h5_read(grp, "measure_G_tau", c.measure_G_tau);
    h5_read(grp, "measure_F_tau", c.measure_F_tau);
    h5::try_read(grp, "measure_G_l", c.measure_G_l);
    h5::try_read(grp, "measure_F_l", c.measure_F_l);
    h5_read(grp, "measure_densities", c.measure_densities);
    h5_read(grp, "measure_average_sign", c.measure_average_sign);
    h5_read(grp, "measure_nn_static", c.measure_nn_static);
    h5_read(grp, "measure_nn_tau", c.measure_nn_tau);
    h5_read(grp, "measure_nn_nu_dlr", c.measure_nn_nu_dlr);
    h5_read(grp, "measure_Sperp_tau", c.measure_Sperp_tau);
    h5_read(grp, "measure_state_hist", c.measure_state_hist);
    h5_read(grp, "measure_g2w", c.measure_g2w);
    h5_read(grp, "measure_g3w", c.measure_g3w);
    h5::try_read(grp, "imag_threshold", c.imag_threshold);
    h5_read(grp, "det_init_size", c.det_init_size);
    h5_read(grp, "det_n_operations_before_check", c.det_n_operations_before_check);
    h5_read(grp, "det_precision_warning", c.det_precision_warning);
    h5_read(grp, "det_precision_error", c.det_precision_error);
    h5_read(grp, "det_singular_threshold", c.det_singular_threshold);
    h5_read(grp, "histogram_max_order", c.histogram_max_order);
    h5_read(grp, "visualize_config", c.visualize_config);
  }

} // namespace triqs_ctseg
