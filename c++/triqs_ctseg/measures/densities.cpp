// Copyright (c) 2022--present, The Simons Foundation
// Copyright (c) 2022--present, Max Planck Institute for Polymer Research, Mainz, Germany
// This file is part of TRIQS/ctseg and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#include "densities.hpp"
#include <itertools/itertools.hpp>
#include <fmt/ostream.h>
#include "../logs.hpp"

namespace triqs_ctseg::measures {

  densities::densities(params_t const &params, work_data_t const &wdata, configuration_t const &config, results_t &results)
     : wdata{wdata}, config{config}, results{results}, measure_densities_{params.measure_densities} {

    n = nda::zeros<double>(config.n_color());

    // Log-binning for auto-correlation: [0] = perturbation order, [1..n_color] = sign * density per color
    log_accs_.reserve(1 + config.n_color());
    for (int i = 0; i < 1 + config.n_color(); ++i) log_accs_.emplace_back(dcomplex{0.0}, -1);

    // Linear binning for density errors (one accumulator per block)
    if (measure_densities_) {
      for (auto const &[bl_name, bl_size] : wdata.gf_struct) {
        dens_bins_.emplace_back(nda::zeros<dcomplex>(bl_size), 128, 1);
      }
    }
  }

  // -------------------------------------

  void densities::accumulate(double s) {

    // Log-bin the perturbation order
    long pert_order = config.Delta_order();
    if (wdata.has_Jperp) pert_order += config.Jperp_order();
    log_accs_[0] << dcomplex(double(pert_order));

    if (measure_densities_) {
      Z += s;
      ++N_;
    }

    // Measure density per color from segment lengths
    for (long bl = 0; auto const &[bl_name, bl_size] : wdata.gf_struct) {
      nda::array<dcomplex, 1> step;
      if (measure_densities_) step.resize(bl_size);

      for (long a = 0; a < bl_size; ++a) {
        long c     = wdata.block_to_color(bl, a);
        double sum = 0;
        for (auto const &seg : config.seglists[c]) sum += double(seg.length());
        n[c] += s * sum;

        auto val = dcomplex(s * sum / tau_t::beta());
        log_accs_[1 + c] << val;
        if (measure_densities_) step(a) = val;
      }

      if (measure_densities_) dens_bins_[bl] << step;
      ++bl;
    }
  }

  // -------------------------------------

  void densities::collect_results(mpi::communicator const &c) {
    using triqs::stat::log_binning;

    // Auto-correlation time from log-binning (always active)
    constexpr int min_samples = 32;
    results.auto_corr_time    = 0.0;
    for (auto &log_acc : log_accs_) {
      auto [mean, errs, taus, effs] = log_acc.mean_errors_and_taus(c, min_samples);
      if (!taus.empty()) { results.auto_corr_time = std::max(results.auto_corr_time, std::real(taus.back())); }
      log_acc = log_binning<dcomplex>{dcomplex{0.0}, -1};
    }
    mpi::broadcast(results.auto_corr_time, c, 0);

    if (measure_densities_) {
      Z = mpi::all_reduce(Z, c);
      n = mpi::all_reduce(n, c);
      n /= (Z * tau_t::beta());

      std::map<std::string, nda::array<double, 1>> densities;
      for (long offset = 0; auto [bl_name, bl_size] : wdata.gf_struct) {
        densities[bl_name] = n[range(offset, offset + bl_size)];
        offset += bl_size;
      }
      if (c.rank() == 0) {
        SPDLOG_INFO("Densities:");
        for (auto &[bl, dens] : densities) SPDLOG_INFO("  {}: {}", bl, fmt::streamed(dens));
        SPDLOG_INFO("Auto-correlation time: {}", results.auto_corr_time);
      }

      results.densities = std::move(densities);

      // Compute error bars from linear binning
      N_ = mpi::all_reduce(N_, c);
      auto norm = std::abs(dcomplex(Z) / dcomplex(N_));
      std::map<std::string, nda::array<double, 1>> densities_errors;
      for (long bl = 0; auto const &[bl_name, bl_size] : wdata.gf_struct) {
        auto [m, err, tau] = dens_bins_[bl].mean_error_and_tau(c);
        densities_errors[bl_name] = nda::array<double, 1>(nda::abs(err) / norm);
        ++bl;
      }
      results.densities_errors = std::move(densities_errors);
    } else {
      if (c.rank() == 0) { SPDLOG_INFO("Auto-correlation time: {}", results.auto_corr_time); }
    }
  }

} // namespace triqs_ctseg::measures
