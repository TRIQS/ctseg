#include "params.hpp"

namespace triqs_ctseg {

  void h5_write(h5::group h5group, std::string subgroup_name, constr_params_t const &c) {

    h5::group grp = subgroup_name.empty() ? h5group : h5group.create_group(subgroup_name);

    h5_write(grp, "beta", c.beta);
    h5_write(grp, "gf_struct", c.gf_struct);
    h5_write(grp, "n_tau", c.n_tau);
    h5_write(grp, "n_tau_bosonic", c.n_tau_bosonic);
  }

  //------------------------------------

  void h5_read(h5::group h5group, std::string subgroup_name, constr_params_t &c) {

    h5::group grp = subgroup_name.empty() ? h5group : h5group.open_group(subgroup_name);

    h5_read(grp, "beta", c.beta);
    h5_read(grp, "gf_struct", c.gf_struct);
    h5_read(grp, "n_tau", c.n_tau);
    h5_read(grp, "n_tau_bosonic", c.n_tau_bosonic);
  }

  //------------------------------------

  void h5_write(h5::group h5group, std::string subgroup_name, solve_params_t const &c) {

    h5::group grp = subgroup_name.empty() ? h5group : h5group.create_group(subgroup_name);

    h5_write(grp, "h_int", c.h_int);
    h5_write(grp, "h_loc0", c.h_loc0);
    h5_write(grp, "n_tau_G", c.n_tau_G);
    h5_write(grp, "n_tau_chi2", c.n_tau_chi2);
    h5_write(grp, "n_cycles", c.n_cycles);
    h5_write(grp, "length_cycle", c.length_cycle);
    h5_write(grp, "n_warmup_cycles", c.n_warmup_cycles);
    h5_write(grp, "random_seed", c.random_seed);
    h5_write(grp, "random_name", c.random_name);
    h5_write(grp, "max_time", c.max_time);
    h5_write(grp, "verbosity", c.verbosity);
    h5_write(grp, "move_insert_segment", c.move_insert_segment);
    h5_write(grp, "move_remove_segment", c.move_remove_segment);
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
    h5_write(grp, "measure_densities", c.measure_densities);
    h5_write(grp, "measure_sign", c.measure_sign);
    h5_write(grp, "measure_nn_static", c.measure_nn_static);
    h5_write(grp, "measure_nn_tau", c.measure_nn_tau);
    h5_write(grp, "measure_sperp_tau", c.measure_sperp_tau);
    h5_write(grp, "measure_state_hist", c.measure_state_hist);
    h5_write(grp, "det_init_size", c.det_init_size);
    h5_write(grp, "det_n_operations_before_check", c.det_n_operations_before_check);
    h5_write(grp, "det_precision_warning", c.det_precision_warning);
    h5_write(grp, "det_precision_error", c.det_precision_error);
    h5_write(grp, "det_singular_threshold", c.det_singular_threshold);
    h5_write(grp, "histogram_max_order", c.histogram_max_order);
  }

  //------------------------------------

  void h5_read(h5::group h5group, std::string subgroup_name, solve_params_t &c) {

    h5::group grp = subgroup_name.empty() ? h5group : h5group.open_group(subgroup_name);

    h5_read(grp, "h_int", c.h_int);
    h5_read(grp, "h_loc0", c.h_loc0);
    h5_read(grp, "n_tau_G", c.n_tau_G);
    h5_read(grp, "n_tau_chi2", c.n_tau_chi2);
    h5_read(grp, "n_cycles", c.n_cycles);
    h5_read(grp, "length_cycle", c.length_cycle);
    h5_read(grp, "n_warmup_cycles", c.n_warmup_cycles);
    h5_read(grp, "random_seed", c.random_seed);
    h5_read(grp, "random_name", c.random_name);
    h5_read(grp, "max_time", c.max_time);
    h5_read(grp, "verbosity", c.verbosity);
    h5_read(grp, "move_insert_segment", c.move_insert_segment);
    h5_read(grp, "move_remove_segment", c.move_remove_segment);
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
    h5_read(grp, "measure_densities", c.measure_densities);
    h5_read(grp, "measure_sign", c.measure_sign);
    h5_read(grp, "measure_nn_static", c.measure_nn_static);
    h5_read(grp, "measure_nn_tau", c.measure_nn_tau);
    h5_read(grp, "measure_sperp_tau", c.measure_sperp_tau);
    h5_read(grp, "measure_state_hist", c.measure_state_hist);
    h5_read(grp, "det_init_size", c.det_init_size);
    h5_read(grp, "det_n_operations_before_check", c.det_n_operations_before_check);
    h5_read(grp, "det_precision_warning", c.det_precision_warning);
    h5_read(grp, "det_precision_error", c.det_precision_error);
    h5_read(grp, "det_singular_threshold", c.det_singular_threshold);
    h5_read(grp, "histogram_max_order", c.histogram_max_order);
  }

} // namespace triqs_ctseg
