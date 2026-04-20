// Copyright (c) 2022--present, The Simons Foundation
// Copyright (c) 2022--present, Max Planck Institute for Polymer Research, Mainz, Germany
// This file is part of TRIQS/ctseg and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#pragma once
#include <optional>
#include <triqs/stat/histograms.hpp>
#include <triqs/gfs.hpp>

using namespace triqs::gfs;

namespace triqs_ctseg {

  // Gather all the results of the CTQMC
  struct results_t {

    /// Single-particle Green's function :math:`G(\tau)`.
    block_gf<imtime> G_tau;

    /// Self-energy improved estimator :math:`F(\tau)`.
    std::optional<block_gf<imtime>> F_tau;

    /// Density-density time correlation function :math:`\langle n_a(\tau) n_b(0) \rangle`.
    std::optional<block2_gf<imtime>> nn_tau;

    /// Density-density frequency correlation function :math:`\langle n_a(i\nu) n_b(-i\nu) \rangle`.
    std::optional<block2_gf<dlr_imfreq>> nn_nu_dlr;

    /// Perpendicular spin-spin correlation function :math:`\langle S_x(\tau) S_x(0) \rangle`.
    std::optional<gf<imtime>> Sperp_tau;

    /// Density-density static correlation function :math:`\langle n_a(0) n_b(0) \rangle`.
    std::optional<std::map<std::pair<std::string, std::string>, nda::matrix<double>>> nn_static;

    /// Density per color, organized by blocks.
    std::optional<std::map<std::string, nda::array<double, 1>>> densities;

    /// Delta perturbation order histogram
    std::optional<std::vector<double>> pert_order_Delta;

    /// Average Delta perturbation order
    std::optional<double> average_order_Delta;

    /// Error bar for average Delta perturbation order
    std::optional<double> average_order_Delta_error;

    /// Jperp perturbation order histogram
    std::optional<std::vector<double>> pert_order_Jperp;

    /// Average Jperp perturbation order
    std::optional<double> average_order_Jperp;

    /// Error bar for average Jperp perturbation order
    std::optional<double> average_order_Jperp_error;

    /// State histogram
    std::optional<nda::vector<double>> state_hist;

    /// Three-point correlation function
    std::optional<block2_gf<prod<imfreq, imfreq>, tensor_valued<4>>> g2w;

    /// Four-point correlation function
    std::optional<block2_gf<prod<imfreq, imfreq, imfreq>, tensor_valued<4>>> g3w;

    /// Auto-correlation time
    double auto_corr_time = 0.0;

    /// Number of warmup cycles actually performed
    int64_t warmup_cycles_done = 0;

    /// The length_cycle value used during accumulation (after auto-determination)
    int length_cycle_used = 0;

    /// Error bars for densities, organized by blocks.
    std::optional<std::map<std::string, nda::array<double, 1>>> densities_errors;

    /// Error bars for density-density static correlations, organized by block pairs.
    std::optional<std::map<std::pair<std::string, std::string>, nda::matrix<double>>> nn_static_errors;

    /// Error bars for state histogram.
    std::optional<nda::vector<double>> state_hist_errors;

    /// Average sign
    double average_sign;

    /// Error bar for average sign.
    std::optional<double> average_sign_error;
  };

  /// writes all containers to hdf5 file
  void h5_write(h5::group h5group, std::string subgroup_name, results_t const &c);

  /// reads all containers from hdf5 file
  void h5_read(h5::group h5group, std::string subgroup_name, results_t &c);

} // namespace triqs_ctseg
