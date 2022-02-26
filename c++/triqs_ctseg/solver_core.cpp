/*******************************************************************************
 * CTSEG: TRIQS hybridization-expansion segment solver
 *
 * Copyright (C) 2013-2018 by T. Ayral, H. Hafermann, P. Delange, M. Ferrero, O.
 *Parcollet Copyright (C) 2019 Simons Foundation author: N. Wentzell
 *
 * CTSEG is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * CTSEG is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * CTSEG. If not, see <http://www.gnu.org/licenses/>.
 ******************************************************************************/
#include "solver_core.hpp"

#include <itertools/itertools.hpp>
#include <triqs/mc_tools/mc_generic.hpp>
#include <triqs/utility/callbacks.hpp>

#include "work_data.hpp"
#include "configuration.hpp"
//#include "measures.hpp"
//#include "moves.hpp"

namespace triqs_ctseg {

solver_core::solver_core(constr_params_t const &p) : constr_params(p) {

  beta = p.beta;

   inputs.delta = block_gf<imtime>(triqs::mesh::imtime{beta, Fermion, p.n_tau}, p.gf_struct);
  inputs.d0t = gf<imtime>({beta, Boson, p.n_tau_jperp}, {1, 1});
  inputs.jperpt = gf<imtime>({beta, Boson, p.n_tau_jperp}, {1, 1});
  
  inputs.delta() = 0;
  inputs.d0t = 0;
  inputs.jperpt = 0;
}

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

  work_data_t wdata{p, inputs};
  configuration_t config{wdata.n_color};

  // ................   QMC  ...................

  auto CTQMC = triqs::mc_tools::mc_generic<double>(p.random_name, p.random_seed,
                                                   p.verbosity);
  
#if 0
  // initialize moves
  if (p.move_insert_segment)
    CTQMC.add_move(move_insert_segment(&params, &config, CTQMC.get_rng()),
                   "segment insertion");
  if (p.move_remove_segment)
    CTQMC.add_move(move_remove_segment(&params, &config, CTQMC.get_rng()),
                   "segment removal");
  if (p.move_move)
    CTQMC.add_move(move_move_segment(&params, &config, CTQMC.get_rng()),
                   "segment move to another line");
  if (p.move_swap_empty_lines)
    CTQMC.add_move(move_swap_empty_lines(&params, &config, CTQMC.get_rng()),
                   "swap empty lines");
  if (jperp_interactions)
    CTQMC.add_move(move_insert_spin_segment(&params, &config, CTQMC.get_rng()),
                   "insert composite segment");
  if (jperp_interactions)
    CTQMC.add_move(move_remove_spin_segment(&params, &config, CTQMC.get_rng()),
                   "remove composite segment");
  if (jperp_interactions)
    CTQMC.add_move(move_swap_bosonic_lines(&params, &config, CTQMC.get_rng()),
                   "swap bosonic lines");
  if (p.move_group_into_spin_segment)
    CTQMC.add_move(
        move_group_into_spin_segment(&params, &config, CTQMC.get_rng()),
        "group into composite segment");
  if (p.move_split_spin_segment)
    CTQMC.add_move(move_split_spin_segment(&params, &config, CTQMC.get_rng()),
                   "split composite segment");
  if (p.move_group_into_spin_segment2)
    CTQMC.add_move(
        move_group_into_spin_segment2(&params, &config, CTQMC.get_rng()),
        "group into composite segment2");
  if (p.move_split_spin_segment2)
    CTQMC.add_move(move_split_spin_segment2(&params, &config, CTQMC.get_rng()),
                   "split composite segment2");

   // FIXME ? REMOVE ? 

  // initialize precomputations shared between measurements
  auto aux_measure_f =
      std::make_shared<precompute_fprefactor>(&params, &config);
  if (p.measure_fw || p.measure_ft || p.measure_fl || p.measure_f2w ||
      p.measure_f3w || p.measure_g2w || p.measure_g3w)
    CTQMC.add_measure_aux(
        aux_measure_f); // MUST BE ADDED BEFORE MEASURE_M SINCE IT IS USED BY IT

  auto aux_measure_M = std::make_shared<precompute_Mw>(
      &params, &p, &config, aux_measure_f, p.n_w_f_vertex, p.n_w_b_vertex);
  if (p.measure_g2w || p.measure_g3w)
    CTQMC.add_measure_aux(aux_measure_M);
  auto aux_measure_n = std::make_shared<precompute_nw>(&params, &config);
  if (p.measure_g2w || p.measure_nnw)
    CTQMC.add_measure_aux(aux_measure_n);

  // initialize measurements
  if (p.measure_gt)
    CTQMC.add_measure(measure_gt(&params, &config, aux_measure_f, gt, ft),
                      "G(tau) measurement");
  if (p.measure_nnt)
    CTQMC.add_measure(measure_nnt(&params, &config, nnt),
                      "nn(tau) measurement");
  if (p.measure_nn)
    CTQMC.add_measure(measure_nn(&params, &config, &nn_matrix),
                      "density measurement");
  if (p.measure_sign)
    CTQMC.add_measure(measure_sign(&average_sign_),
                      "measurement of the MC sign");
 
   // FIXME  ? Keep ? ?
 
  if (p.measure_hist)
    CTQMC.add_measure(measure_hist(&config, &hist), "histogram measurement");
 
  if (p.measure_hist_composite)
    CTQMC.add_measure(measure_hist_composite(&config, &hist_composite),
                      "histogram of composite order measurement");
  if (p.measure_statehist)
    CTQMC.add_measure(measure_statehist(&params, &config, &state_hist),
                      "impurity state histogram measurement");

#endif

  // Run and collect results
  auto _solve_status =
      CTQMC.warmup_and_accumulate(p.n_warmup_cycles, p.n_cycles, p.length_cycle,
                                  triqs::utility::clock_callback(p.max_time));
  CTQMC.collect_results(c);

} // solve

} // namespace triqs_ctseg
