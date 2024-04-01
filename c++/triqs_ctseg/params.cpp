#include "params.hpp"

void h5_write(h5::group h5group, std::string subgroup_name, constr_params_t const &c) {

  h5::group grp = subgroup_name.empty() ? h5group : h5group.create_group(subgroup_name);

  h5_write(grp, "beta", c.beta);
  h5_write(grp, "gf_struct", c.gf_struct);
  h5_write(grp, "n_tau", c.n_tau);
  h5_write(grp, "n_tau_k", c.n_tau_k);
}

//------------------------------------

void h5_read(h5::group h5group, std::string subgroup_name, constr_params_t &c) {

  h5::group grp = subgroup_name.empty() ? h5group : h5group.open_group(subgroup_name);

  h5_read(grp, "beta", c.beta);
  h5_read(grp, "gf_struct", c.gf_struct);
  h5_read(grp, "n_tau", c.n_tau);
  h5_read(grp, "n_tau_k", c.n_tau_k);
}

//------------------------------------

void h5_write(h5::group h5group, std::string subgroup_name, solve_params_t const &c) {

  h5::group grp = subgroup_name.empty() ? h5group : h5group.create_group(subgroup_name);

  h5_write(grp, "h_int", c.h_int);
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
  h5_write(grp, "measure_perturbation_order_histograms", c.measure_perturbation_order_histograms);
  h5_write(grp, "measure_gt", c.measure_gt);
  h5_write(grp, "measure_ft", c.measure_ft);
  h5_write(grp, "measure_n", c.measure_n);
  h5_write(grp, "measure_sign", c.measure_sign);
  h5_write(grp, "measure_nn", c.measure_nn);
  h5_write(grp, "measure_nnt", c.measure_nnt);
  h5_write(grp, "measure_sperpt", c.measure_sperpt);
  h5_write(grp, "measure_statehist", c.measure_statehist);
  h5_write(grp, "hartree_shift", c.hartree_shift);
  h5_write(grp, "det_init_size", c.det_init_size);
  h5_write(grp, "det_n_operations_before_check", c.det_n_operations_before_check);
  h5_write(grp, "det_precision_warning", c.det_precision_warning);
  h5_write(grp, "det_precision_error", c.det_precision_error);
  h5_write(grp, "det_singular_threshold", c.det_singular_threshold);
  
}

//------------------------------------

void h5_read(h5::group h5group, std::string subgroup_name, solve_params_t &c) {

  h5::group grp = subgroup_name.empty() ? h5group : h5group.open_group(subgroup_name);

  h5_read(grp, "h_int", c.h_int);
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
  h5_read(grp, "measure_perturbation_order_histograms", c.measure_perturbation_order_histograms);
  h5_read(grp, "measure_gt", c.measure_gt);
  h5_read(grp, "measure_ft", c.measure_ft);
  h5_read(grp, "measure_n", c.measure_n);
  h5_read(grp, "measure_sign", c.measure_sign);
  h5_read(grp, "measure_nn", c.measure_nn);
  h5_read(grp, "measure_nnt", c.measure_nnt);
  h5_read(grp, "measure_sperpt", c.measure_sperpt);
  h5_read(grp, "measure_statehist", c.measure_statehist);
  h5_read(grp, "hartree_shift", c.hartree_shift);
  h5_read(grp, "det_init_size", c.det_init_size);
  h5_read(grp, "det_n_operations_before_check", c.det_n_operations_before_check);
  h5_read(grp, "det_precision_warning", c.det_precision_warning);
  h5_read(grp, "det_precision_error", c.det_precision_error);
  h5_read(grp, "det_singular_threshold", c.det_singular_threshold);

}

