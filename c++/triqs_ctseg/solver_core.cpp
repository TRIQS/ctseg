#include "solver_core.hpp"

#include <triqs/mc_tools/mc_generic.hpp>
#include <triqs/utility/callbacks.hpp>

#include "work_data.hpp"
#include "configuration.hpp"

// FIXME : remove those includers ??
#include "measures.hpp"
#include "moves.hpp"

#include "logs.hpp"

solver_core::solver_core(constr_params_t const &p) : constr_params(p) {

  beta = p.beta;

  inputs.delta  = block_gf<imtime>(triqs::mesh::imtime{beta, Fermion, p.n_tau}, p.gf_struct);
  inputs.d0t    = gf<imtime>({beta, Boson, p.n_tau_k}, {1, 1});
  inputs.jperpt = gf<imtime>({beta, Boson, p.n_tau_jperp}, {1, 1});

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
                 " ║ ╠╦╝║║═╬╗╚═╗  │   │ └─┐├┤ │ ┬\n"
                 " ╩ ╩╚═╩╚═╝╚╚═╝  └─┘ ┴ └─┘└─┘└─┘\n";

  // parameters
  last_solve_params = solve_params;
  // FIXME ? keep it ?
  // Merge constr_params and solve_params
  params_t p(constr_params, solve_params);

  // ................   wdata & config  ...................

  work_data_t wdata{p, inputs, c};
  configuration_t config{wdata.n_color};
  results.K_tau = wdata.K;

  // ................   QMC  ...................

  auto CTQMC = triqs::mc_tools::mc_generic<double>(p.random_name, p.random_seed, p.verbosity);

  // initialize moves
  if (p.move_insert_segment) CTQMC.add_move(moves::insert_segment{wdata, config, CTQMC.get_rng()}, "insert");
  if (p.move_remove_segment) CTQMC.add_move(moves::remove_segment{wdata, config, CTQMC.get_rng()}, "remove");
  if (p.move_move_segment) CTQMC.add_move(moves::move_segment{wdata, config, CTQMC.get_rng()}, "move");
  if (p.move_split_segment) CTQMC.add_move(moves::split_segment{wdata, config, CTQMC.get_rng()}, "split");
  if (p.move_regroup_segment) CTQMC.add_move(moves::regroup_segment{wdata, config, CTQMC.get_rng()}, "regroup");
  if (p.move_insert_spin_segment)
    CTQMC.add_move(moves::insert_spin_segment{wdata, config, CTQMC.get_rng()}, "spin insert");
  if (p.move_remove_spin_segment)
    CTQMC.add_move(moves::remove_spin_segment{wdata, config, CTQMC.get_rng()}, "spin remove");
  if (p.move_swap_spin_lines) CTQMC.add_move(moves::swap_spin_lines{wdata, config, CTQMC.get_rng()}, "spin swap");

  // unused moves - for testing purposes
  /* if (p.move_insert_segment_v2) CTQMC.add_move(moves::insert_segment_v2{wdata, config, CTQMC.get_rng()}, "insert v2");
  if (p.move_remove_segment_v2) CTQMC.add_move(moves::remove_segment_v2{wdata, config, CTQMC.get_rng()}, "remove v2");
  if (p.move_split_segment_v2) CTQMC.add_move(moves::split_segment_v2{wdata, config, CTQMC.get_rng()}, "split v2");
  if (p.move_regroup_segment_v2)
    CTQMC.add_move(moves::regroup_segment_v2{wdata, config, CTQMC.get_rng()}, "regroup v2"); */

  // initialize measurements
  if (p.measure_gt) CTQMC.add_measure(measures::g_f_tau{p, wdata, config, results}, "G(tau)");
  if (p.measure_n) CTQMC.add_measure(measures::density{p, wdata, config, results}, "Density");
  if (p.measure_nnt) CTQMC.add_measure(measures::nn_tau{p, wdata, config, results}, "nn(tau)");

  if (p.measure_perturbation_order_histograms)
    CTQMC.add_measure(measures::perturbation_order_histo{p, wdata, config, results}, "Perturbation orders");

  // Run and collect results
  auto _solve_status = CTQMC.warmup_and_accumulate(p.n_warmup_cycles, p.n_cycles, p.length_cycle,
                                                   triqs::utility::clock_callback(p.max_time));
  CTQMC.collect_results(c);
  //SPDLOG_INFO("Final config {}", config);

}; // solve
