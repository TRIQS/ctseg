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

  /// Container for all results accumulated by the CTQMC simulation.
  struct results_t {

    /// Single-particle Green's function \f$ G(\tau) \f$.
    block_gf<imtime> G_tau;

    /// Self-energy improved estimator \f$ F(\tau) \f$.
    std::optional<block_gf<imtime>> F_tau;

    /// Density-density time correlation function \f$ \langle n_a(\tau) n_b(0) \rangle \f$.
    std::optional<block2_gf<imtime>> nn_tau;

    /// Density-density frequency correlation function \f$ \langle n_a(i\nu) n_b(-i\nu) \rangle \f$.
    std::optional<block2_gf<dlr_imfreq>> nn_nu_dlr;

    /// Perpendicular spin-spin correlation function \f$ \langle S_x(\tau) S_x(0) \rangle \f$.
    std::optional<gf<imtime>> Sperp_tau;

    /// Oriented transverse spin correlation \f$ \langle S^-(\tau) S^+(0) \rangle \f$.
    std::optional<gf<imtime>> Sminus_Splus_tau;

    /// Oriented transverse spin correlation \f$ \langle S^+(\tau) S^-(0) \rangle \f$.
    std::optional<gf<imtime>> Splus_Sminus_tau;

    /// Density-density static correlation function \f$ \langle n_a(0) n_b(0) \rangle \f$.
    std::optional<std::map<std::pair<std::string, std::string>, nda::matrix<double>>> nn_static;

    /// Density per color, organized by blocks.
    std::optional<std::map<std::string, nda::array<double, 1>>> densities;

    /// Delta perturbation order histogram.
    std::optional<std::vector<double>> pert_order_Delta;

    /// Average Delta perturbation order.
    std::optional<double> average_order_Delta;

    /// Jperp perturbation order histogram.
    std::optional<std::vector<double>> pert_order_Jperp;

    /// Average Jperp perturbation order.
    std::optional<double> average_order_Jperp;

    /// State histogram.
    std::optional<nda::vector<double>> state_hist;

    /// Retarded source-field/density static correlation in color space.
    std::optional<nda::matrix<double>> dyn_phi_n;

    /// Retarded source-field/source-field static correlation in color space.
    std::optional<nda::matrix<double>> dyn_phi_phi;

    /// Three-point correlation function.
    std::optional<block2_gf<prod<imfreq, imfreq>, tensor_valued<4>>> g2w;

    /// Four-point correlation function.
    std::optional<block2_gf<prod<imfreq, imfreq, imfreq>, tensor_valued<4>>> g3w;

    /// Average sign.
    double average_sign;

    static std::string hdf5_format() { return "CTSEG_Results"; }
  };

  /// writes all containers to hdf5 file
  void h5_write(h5::group h5group, std::string subgroup_name, results_t const &c);

  /// reads all containers from hdf5 file
  void h5_read(h5::group h5group, std::string subgroup_name, results_t &c);

} // namespace triqs_ctseg
