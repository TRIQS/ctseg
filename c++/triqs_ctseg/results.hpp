// Copyright (c) 2022--present, The TRIQS/ctseg Authors and their Assignees
// This file is part of TRIQS/ctseg and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE.txt in the root of this distribution for details.

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

    /// Jperp perturbation order histogram
    std::optional<std::vector<double>> pert_order_Jperp;

    /// Average Jperp perturbation order
    std::optional<double> average_order_Jperp;

    /// State histogram
    std::optional<nda::vector<double>> state_hist;

    /// Four-point correlation function
    std::optional<block2_gf<prod<imfreq, imfreq, imfreq>, tensor_valued<4>>> g3w;

    /// Four-point correlation function improved estimator
    std::optional<block2_gf<prod<imfreq, imfreq, imfreq>, tensor_valued<4>>> f3w;

    /// Average sign
    double average_sign;
  };

  /// writes all containers to hdf5 file
  void h5_write(h5::group h5group, std::string subgroup_name, results_t const &c);

  /// reads all containers from hdf5 file
  void h5_read(h5::group h5group, std::string subgroup_name, results_t &c);

} // namespace triqs_ctseg
