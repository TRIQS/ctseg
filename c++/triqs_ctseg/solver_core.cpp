#include <triqs/mc_tools/mc_generic.hpp>
#include <triqs/utility/callbacks.hpp>

#include "solver_core.hpp"
#include "work_data.hpp"
#include "configuration.hpp"
#include "measures.hpp"
#include "moves.hpp"
#include "logs.hpp"

// ---------------------------------------------------------------------------

solver_core::solver_core(constr_params_t const &p) : constr_params(p) {

  beta = p.beta;
  tau_t::set_beta(beta);
  int n_color = count_colors(p.gf_struct);

  inputs.delta  = block_gf<imtime>(triqs::mesh::imtime{beta, Fermion, p.n_tau}, p.gf_struct);
  inputs.d0t    = gf<imtime>({beta, Boson, p.n_tau_k}, {n_color, n_color});
  inputs.jperpt = gf<imtime>({beta, Boson, p.n_tau_k}, {1, 1});

  inputs.delta()  = 0;
  inputs.d0t()    = 0;
  inputs.jperpt() = 0;
};

// ---------------------------------------------------------------------------

void solver_core::solve(solve_params_t const &solve_params) {

  // http://patorjk.com/software/taag/#p=display&f=Calvin%20S&t=TRIQS%20ctint
  if (c.rank() == 0)
    std::cout << "\n"
                 "╔╦╗╦═╗╦╔═╗ ╔═╗  ┌─┐┌┬┐┌─┐┌─┐┌─┐\n"
                 " ║ ╠╦╝║║═╬╗╚═╗  │   │ └─┐├┤ │ ┬    V2-J \n"
                 " ╩ ╩╚═╩╚═╝╚╚═╝  └─┘ ┴ └─┘└─┘└─┘\n";

  // parameters
  last_solve_params = solve_params;
  // Merge constr_params and solve_params
  params_t p(constr_params, solve_params);

  // ................   wdata & config  ...................

  work_data_t wdata{p, inputs, c};
  configuration_t config{wdata.n_color};
  results.K_tau      = wdata.K;
  results.Kprime_tau = wdata.Kprime;

  // ................   QMC  ...................

  auto CTQMC = triqs::mc_tools::mc_generic<double>(p.random_name, p.random_seed, p.verbosity);

  // initialize moves
  if (p.move_insert_segment) CTQMC.add_move(moves::insert_segment{wdata, config, CTQMC.get_rng()}, "insert");
  if (p.move_remove_segment) CTQMC.add_move(moves::remove_segment{wdata, config, CTQMC.get_rng()}, "remove");
  if (p.move_move_segment) CTQMC.add_move(moves::move_segment{wdata, config, CTQMC.get_rng()}, "move");
  if (p.move_split_segment) CTQMC.add_move(moves::split_segment{wdata, config, CTQMC.get_rng()}, "split");
  if (p.move_regroup_segment) CTQMC.add_move(moves::regroup_segment{wdata, config, CTQMC.get_rng()}, "regroup");

  if (wdata.has_jperp) {
    if (p.move_insert_spin_segment)
      CTQMC.add_move(moves::insert_spin_segment{wdata, config, CTQMC.get_rng()}, "spin insert");

    if (p.move_remove_spin_segment)
      CTQMC.add_move(moves::remove_spin_segment{wdata, config, CTQMC.get_rng()}, "spin remove");

    if (p.move_split_spin_segment)
      CTQMC.add_move(moves::split_spin_segment{wdata, config, CTQMC.get_rng()}, "spin split");

    if (p.move_regroup_spin_segment)
      CTQMC.add_move(moves::regroup_spin_segment{wdata, config, CTQMC.get_rng()}, "spin regroup");

    if (p.move_swap_spin_lines) CTQMC.add_move(moves::swap_spin_lines{wdata, config, CTQMC.get_rng()}, "spin swap");
  }

  // initialize measurements
  if (p.measure_gt) CTQMC.add_measure(measures::g_f_tau{p, wdata, config, results}, "G(tau)");
  if (p.measure_n) CTQMC.add_measure(measures::density{p, wdata, config, results}, "Density");
  if (p.measure_nn) CTQMC.add_measure(measures::nn_static{p, wdata, config, results}, "nn(0)");
  if (p.measure_nnt) CTQMC.add_measure(measures::nn_tau{p, wdata, config, results}, "nn(tau)");
  if (p.measure_sperpt) CTQMC.add_measure(measures::sperp_tau{p, wdata, config, results}, "<s_x s_x>(tau)");
  if (p.measure_perturbation_order_histograms)
    CTQMC.add_measure(measures::perturbation_order_histo{p, wdata, config, results}, "Perturbation orders");

  // Run and collect results
  CTQMC.warmup_and_accumulate(p.n_warmup_cycles, p.n_cycles, p.length_cycle,
                              triqs::utility::clock_callback(p.max_time));
  CTQMC.collect_results(c);

} // solve

// -------------- h5 -----------------------

#define STR(x) #x
#define STRINGIZE(x) STR(x)

// Function that writes the solver_core to hdf5 file
void h5_write(h5::group h5group, std::string subgroup_name, solver_core const &s) {
  auto grp = h5group.create_group(subgroup_name);
  h5_write_attribute(grp, "Format", solver_core::hdf5_format());
  h5_write_attribute(grp, "TRIQS_GIT_HASH", std::string(STRINGIZE(TRIQS_GIT_HASH)));
  h5_write_attribute(grp, "CTSEG_GIT_HASH", std::string(STRINGIZE(CTSEG_GIT_HASH)));
  h5_write(grp, "constr_params", s.constr_params);
  // FIXME : implement
  //h5_write(grp, "last_solve_params", s.last_solve_params);
  h5_write(grp, "inputs", s.inputs);
  h5_write(grp, "results", s.results);
}

// Function that read all containers to hdf5 file
solver_core solver_core::h5_read_construct(h5::group h5group, std::string subgroup_name) {
  auto grp           = h5group.open_group(subgroup_name);
  auto constr_params = h5_read<constr_params_t>(grp, "constr_params");
  auto s             = solver_core{constr_params};
  //h5_read(grp, "last_solve_params", s.last_solve_params);
  h5_read(grp, "inputs", s.inputs);
  h5_read(grp, "results", s.results);
  return s;
}
