#include "solver_core.hpp"

#include <triqs/mc_tools/mc_generic.hpp>
#include <triqs/utility/callbacks.hpp>

#include "work_data.hpp"
#include "configuration.hpp"
#include "measures.hpp"
#include "moves.hpp"

using namespace moves; 
using namespace measures; 

solver_core::solver_core(constr_params_t const &p) : constr_params(p) {

  beta = p.beta;

  inputs.delta  = block_gf<imtime>(triqs::mesh::imtime{beta, Fermion, p.n_tau}, p.gf_struct);
  inputs.d0t    = gf<imtime>({beta, Boson, p.n_tau_jperp}, {1, 1});
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

  // ................   QMC  ...................

  auto CTQMC = triqs::mc_tools::mc_generic<double>(p.random_name, p.random_seed, p.verbosity);

  // initialize moves
  if (p.move_insert_segment)
    CTQMC.add_move(insert_segment{wdata, config, CTQMC.get_rng()},
                   "segment insertion");
  if (p.move_remove_segment)
    CTQMC.add_move(remove_segment(wdata, config, CTQMC.get_rng()),
                   "segment removal");
  if (p.move_move)
    CTQMC.add_move(move_segment(wdata, config, CTQMC.get_rng()),
                   "segment move to another line");

  if (p.move_split_segment)
    CTQMC.add_move(split_segment(wdata, config, CTQMC.get_rng()),
                   "split a segment");

  if (p.move_regroup_segment)
    CTQMC.add_move(regroup_segment(wdata, config, CTQMC.get_rng()),
                   "regroup two segments");                                     

  // initialize measurements
  if (p.measure_gt)
    CTQMC.add_measure(g_f_tau(p, wdata, config, results),
                      "G(tau) measurement");
  if (p.measure_nnt)
    CTQMC.add_measure(nn_tau(p, wdata, config, results),
                      "nn(tau) measurement");
  //if (p.measure_nn)
  //  CTQMC.add_measure(measure_nn(&params, &config, &nn_matrix),
  //                   "density measurement");
 
   // FIXME  ? Keep ? ?
  #if 0
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
  auto _solve_status = CTQMC.warmup_and_accumulate(p.n_warmup_cycles, p.n_cycles, p.length_cycle,
                                                   triqs::utility::clock_callback(p.max_time));
  CTQMC.collect_results(c);

}; // solve
