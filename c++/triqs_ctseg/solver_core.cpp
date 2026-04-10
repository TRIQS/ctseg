// Copyright (c) 2022--present, The Simons Foundation
// Copyright (c) 2022--present, Max Planck Institute for Polymer Research, Mainz, Germany
// This file is part of TRIQS/ctseg and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#include <triqs/mc_tools/mc_generic.hpp>
#include <triqs/utility/callbacks.hpp>
#include <triqs/stat/log_binning.hpp>
#include <chrono>

#include "solver_core.hpp"
#include "work_data.hpp"
#include "configuration.hpp"
#include "measures.hpp"
#include "moves.hpp"
#include "logs.hpp"

namespace triqs_ctseg {

  // ---------------------------------------------------------------------------

  solver_core::solver_core(constr_params_t const &p) : constr_params(p) {

    beta = p.beta;
    tau_t::set_beta(beta);

    inputs.Delta  = block_gf<imtime>({beta, Fermion, p.n_tau}, p.gf_struct);
    inputs.D0t    = make_block2_gf<imtime>({beta, Boson, p.n_tau_bosonic}, p.gf_struct);
    inputs.Jperpt = gf<imtime>({beta, Boson, p.n_tau_bosonic}, {1, 1});

    inputs.Delta()  = 0;
    inputs.D0t()    = 0;
    inputs.Jperpt() = 0;
  };

  // ---------------------------------------------------------------------------

  void solver_core::solve(solve_params_t const &solve_params_input) {

    // http://patorjk.com/software/taag/#p=display&f=Calvin%20S&t=TRIQS%20ctseg
    if (c.rank() == 0)
      std::cout << "\n"
                   "╔╦╗╦═╗╦╔═╗ ╔═╗  ┌─┐┌┬┐┌─┐┌─┐┌─┐\n"
                   " ║ ╠╦╝║║═╬╗╚═╗  │   │ └─┐├┤ │ ┬\n"
                   " ╩ ╩╚═╩╚═╝╚╚═╝  └─┘ ┴ └─┘└─┘└─┘\n";

    // ................ Parameters .................
    // Store the solve_params
    solve_params = solve_params_input;
    // Set tau mesh parameters for results to default if not supplied
    if (solve_params_input.n_tau_G == 0) solve_params.n_tau_G = constr_params.n_tau;
    if (solve_params_input.n_tau_G == 0) solve_params.n_tau_chi2 = constr_params.n_tau_bosonic;
    // Merge constr_params and solve_params
    params_t p(constr_params, solve_params);

    // ................   Work data & Configuration  ...................

    // Initialize work data
    work_data_t wdata{p, inputs, c};
    // Initialize configuration
    configuration_t config{wdata.n_color};
    // Start from a non-empty configuration when Delta(tau) = 0
    if (not wdata.has_Delta) { config.seglists[0].push_back(segment_t::full_line()); }

    // ................   QMC  ...................

    auto CTQMC = triqs::mc_tools::mc_generic<double>(p.random_name, p.random_seed, p.verbosity);

    // Initialize moves
    if (wdata.has_Delta) {
      if (p.move_insert_segment) CTQMC.add_move(moves::insert_segment{wdata, config, CTQMC.get_rng()}, "insert");
      if (p.move_remove_segment) CTQMC.add_move(moves::remove_segment{wdata, config, CTQMC.get_rng()}, "remove");
      if (p.move_double_insert_segment) CTQMC.add_move(moves::double_insert_segment{wdata, config, CTQMC.get_rng()}, "double insert");
      if (p.move_double_remove_segment) CTQMC.add_move(moves::double_remove_segment{wdata, config, CTQMC.get_rng()}, "double remove");
      if (p.move_move_segment) CTQMC.add_move(moves::move_segment{wdata, config, CTQMC.get_rng()}, "move");
      if (p.move_split_segment) CTQMC.add_move(moves::split_segment{wdata, config, CTQMC.get_rng()}, "split");
      if (p.move_regroup_segment) CTQMC.add_move(moves::regroup_segment{wdata, config, CTQMC.get_rng()}, "regroup");
    }

    if (wdata.has_Jperp) {
      if (p.move_insert_spin_segment)
        CTQMC.add_move(moves::insert_spin_segment{wdata, config, CTQMC.get_rng()}, "spin insert");

      if (p.move_remove_spin_segment)
        CTQMC.add_move(moves::remove_spin_segment{wdata, config, CTQMC.get_rng()}, "spin remove");
    }

    if (wdata.has_Jperp and wdata.has_Delta) {
      if (p.move_split_spin_segment)
        CTQMC.add_move(moves::split_spin_segment{wdata, config, CTQMC.get_rng()}, "spin split");

      if (p.move_regroup_spin_segment)
        CTQMC.add_move(moves::regroup_spin_segment{wdata, config, CTQMC.get_rng()}, "spin regroup");
    }

    if (wdata.has_Jperp) {
      if (p.move_swap_spin_lines) CTQMC.add_move(moves::swap_spin_lines{wdata, config, CTQMC.get_rng()}, "spin swap");
    }

    // ========== Phase 1: Warmup ==========

    int warmup_cycle_length = (p.length_cycle >= 0) ? p.length_cycle : 100;
    bool auto_warmup        = (p.n_warmup_cycles < 0);

    if (auto_warmup) {
      if (c.rank() == 0) spdlog::info("Warming up (auto) ...");
      CTQMC.set_verbosity(0);

      // Automatic warmup: run until perturbation order stabilizes
      double mean_k      = 0.0;
      double prev_mean_k = 0.0;
      double sign_sum    = 0.0;
      int64_t n_acc      = 0;
      int n_stable       = 0;
      bool converged     = false;
      double next_print  = 2.0;

      constexpr int check_interval    = 100;
      constexpr int min_warmup        = 100;
      constexpr double rtol           = 0.03;
      constexpr int n_stable_required = 3;

      auto clock_cb = triqs::utility::clock_callback(p.max_time);
      auto t0       = std::chrono::steady_clock::now();

      auto after_duty = [&]() {
        long pert_order = config.Delta_order();
        if (wdata.has_Jperp) pert_order += config.Jperp_order();
        double k = static_cast<double>(pert_order);
        ++n_acc;
        mean_k += (k - mean_k) / n_acc; // online mean
        sign_sum += CTQMC.get_sign();
      };

      auto stop_cb = [&]() -> bool {
        if (clock_cb()) return true;
        if (n_acc < min_warmup || n_acc % check_interval != 0) return false;

        // Periodic status print
        double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        if (elapsed > next_print) {
          next_print = 1.25 * elapsed + 2.0;
          if (c.rank() == 0)
            spdlog::info("  mean k = {}, sign = {} ({} cycles)", mean_k, sign_sum / static_cast<double>(n_acc), n_acc);
        }

        // Compute relative change of running mean
        double rel_change = std::abs(mean_k - prev_mean_k) / std::max(std::abs(mean_k), 1.0);
        prev_mean_k       = mean_k;

        if (rel_change < rtol)
          ++n_stable;
        else
          n_stable = 0;

        bool local_converged = (n_stable >= n_stable_required);
        int global_converged = mpi::all_reduce(static_cast<int>(local_converged), c, MPI_MIN);
        converged            = (global_converged == 1);
        return converged;
      };

      typename decltype(CTQMC)::run_param_t rp;
      rp.ncycles         = p.max_warmup_cycles;
      rp.cycle_length    = warmup_cycle_length;
      rp.stop_callback   = stop_cb;
      rp.after_cycle_duty = after_duty;
      rp.comm            = c;
      rp.enable_measures = false;
      CTQMC.run(rp);

      if (!converged) {
        if (c.rank() == 0) spdlog::warn("Warmup did not converge after {} cycles", p.max_warmup_cycles);
      } else {
        if (c.rank() == 0) spdlog::info("  mean k = {} -> converged", mean_k);
      }

      CTQMC.set_verbosity(p.verbosity);
    } else {
      if (c.rank() == 0) spdlog::info("Warming up ...");
      CTQMC.run(p.n_warmup_cycles, warmup_cycle_length, triqs::utility::clock_callback(p.max_time), /* enable_measures */ false, c);
    }
    results.warmup_cycles_done = CTQMC.get_current_cycle_number();

    // ========== Phase 2: length_cycle calibration ==========

    int effective_length_cycle = p.length_cycle;
    bool auto_length_cycle     = (p.length_cycle < 0);
    if (auto_length_cycle) {
      if (c.rank() == 0) spdlog::info("Calibrating length_cycle ...");
      CTQMC.set_verbosity(0);

      // Register densities measure for calibration (gives auto_corr_time including density correlations)
      CTQMC.add_measure(measures::densities{p, wdata, config, results}, "calibration densities");

      // Log-binning on perturbation order for convergence detection
      triqs::stat::log_binning<dcomplex> k_acc(dcomplex{0.0}, -1);
      int64_t calib_count = 0;
      double prev_tau     = -1.0;
      int n_stable        = 0;
      double next_print   = 2.0;

      constexpr int check_interval    = 500;
      constexpr int min_calib         = 1000;
      constexpr int max_calib_cycles   = 100000;
      constexpr double tau_rtol       = 0.1;
      constexpr int n_stable_required = 3;

      auto clock_cb = triqs::utility::clock_callback(p.max_time);
      auto t0       = std::chrono::steady_clock::now();

      auto after_duty = [&]() {
        long pert_order = config.Delta_order();
        if (wdata.has_Jperp) pert_order += config.Jperp_order();
        k_acc << dcomplex(double(pert_order));
        ++calib_count;
      };

      auto stop_cb = [&]() -> bool {
        if (clock_cb()) return true;
        if (calib_count < min_calib || calib_count % check_interval != 0) return false;

        auto [mean, errs, taus, effs] = k_acc.mean_errors_and_taus(c);
        double tau = taus.empty() ? 0.0 : std::real(taus.back());

        // Periodic status print
        double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        if (elapsed > next_print) {
          next_print = 1.25 * elapsed + 2.0;
          if (c.rank() == 0) spdlog::info("  tau_ac = {} ({} cycles)", tau, calib_count);
        }

        if (prev_tau >= 0) {
          double rel_change = std::abs(tau - prev_tau) / std::max(tau, 1.0);
          if (rel_change < tau_rtol)
            ++n_stable;
          else
            n_stable = 0;
        }
        prev_tau = tau;
        return (n_stable >= n_stable_required);
      };

      typename decltype(CTQMC)::run_param_t rp;
      rp.ncycles          = max_calib_cycles;
      rp.cycle_length     = 1;
      rp.stop_callback    = stop_cb;
      rp.after_cycle_duty = after_duty;
      rp.comm             = c;
      rp.enable_measures  = true;
      CTQMC.run(rp);

      // Collect results to compute auto_corr_time from the densities measure
      CTQMC.collect_results(c);
      double tau_raw = results.auto_corr_time;

      // Set length_cycle so that effective autocorrelation ~ target_auto_corr_time
      effective_length_cycle = std::max(1, static_cast<int>(std::ceil(tau_raw / p.target_auto_corr_time)));
      effective_length_cycle = std::min(effective_length_cycle, p.max_length_cycle);

      if (c.rank() == 0) spdlog::info("  tau_ac = {} -> length_cycle = {}", tau_raw, effective_length_cycle);

      // Clean up calibration phase
      CTQMC.clear_measures();
      results = results_t{};
      results.warmup_cycles_done = CTQMC.get_current_cycle_number();
      CTQMC.set_verbosity(p.verbosity);
    }
    results.length_cycle_used = effective_length_cycle;

    // ========== Phase 3: Accumulation ==========

    // Initialize measurements
    if (p.measure_G_tau) CTQMC.add_measure(measures::G_F_tau{p, wdata, config, results}, "G(tau)/F(tau)");
    CTQMC.add_measure(measures::densities{p, wdata, config, results}, "Densities");
    if (p.measure_average_sign) CTQMC.add_measure(measures::average_sign{p, wdata, config, results}, "Average Sign");
    if (p.measure_nn_static) CTQMC.add_measure(measures::nn_static{p, wdata, config, results}, "<nn>");
    if (p.measure_nn_tau) CTQMC.add_measure(measures::nn_tau{p, wdata, config, results}, "<n(tau)n(0)>");
    if (p.measure_nn_nu_dlr) CTQMC.add_measure(measures::nn_nu_dlr{p, wdata, config, results}, "<n(nu)n(-nu)>");
    if (p.measure_Sperp_tau) CTQMC.add_measure(measures::Sperp_tau{p, wdata, config, results}, "<S_x(tau)S_x(0)>");
    if (p.measure_pert_order) {
      if (wdata.has_Delta) {
        CTQMC.add_measure(measures::pert_order{[&]() { return config.Delta_order(); }, results.pert_order_Delta,
                                               results.average_order_Delta},
                          "Perturbation order Delta");
      }
      if (wdata.has_Jperp) {
        CTQMC.add_measure(measures::pert_order{[&]() { return config.Jperp_order(); }, results.pert_order_Jperp,
                                               results.average_order_Jperp},
                          "Perturbation order Jperp");
      }
    }
    if (p.measure_state_hist) CTQMC.add_measure(measures::state_hist{p, wdata, config, results}, "State histograms");
    if (p.measure_g2w || p.measure_g3w)
      CTQMC.add_measure(measures::four_point{p, wdata, config, results}, "Four-point correlation function");
    if (p.visualize_config) CTQMC.add_measure(measures::visualize_config{config}, "Visualizing configurations");

    // Run accumulation and collect results
    CTQMC.run(p.n_cycles, effective_length_cycle, triqs::utility::clock_callback(p.max_time), /* enable_measures */ true, c);
    CTQMC.collect_results(c);

    // Report summary
    if (c.rank() == 0) {
      spdlog::info("Average sign: {}", results.average_sign);
      if (results.average_order_Delta)
        spdlog::info("Average perturbation order in Delta: {:.3f}", results.average_order_Delta.value());
      if (results.average_order_Jperp)
        spdlog::info("Average perturbation order in Jperp: {:.3f}", results.average_order_Jperp.value());
      spdlog::info("Auto-correlation time: {}", results.auto_corr_time);
      spdlog::info("Warmup cycles: {}{}", results.warmup_cycles_done, auto_warmup ? " (auto)" : "");
      spdlog::info("Length cycle: {}{}", results.length_cycle_used, auto_length_cycle ? " (auto)" : "");
    }

  } // solve

  // ----------------- Save to h5 file -----------------------

#define STR(x) #x
#define STRINGIZE(x) STR(x)

  // Function that writes the solver_core to hdf5 file
  void h5_write(h5::group h5group, std::string subgroup_name, solver_core const &s) {
    auto grp = h5group.create_group(subgroup_name);
    h5_write_attribute(grp, "Format", solver_core::hdf5_format());
    h5_write_attribute(grp, "TRIQS_GIT_HASH", std::string(STRINGIZE(TRIQS_GIT_HASH)));
    h5_write_attribute(grp, "CTSEG_GIT_HASH", std::string(STRINGIZE(CTSEG_GIT_HASH)));
    h5_write(grp, "constr_params", s.constr_params);
    h5_write(grp, "solve_params", s.solve_params);
    h5_write(grp, "inputs", s.inputs);
    h5_write(grp, "results", s.results);
  }

  // Function that reads all containers in hdf5 file
  solver_core solver_core::h5_read_construct(h5::group h5group, std::string subgroup_name) {
    auto grp           = h5group.open_group(subgroup_name);
    auto constr_params = h5_read<constr_params_t>(grp, "constr_params");
    auto s             = solver_core{constr_params};
    h5_read(grp, "solve_params", s.solve_params);
    h5_read(grp, "inputs", s.inputs);
    h5_read(grp, "results", s.results);
    return s;
  }

} // namespace triqs_ctseg
