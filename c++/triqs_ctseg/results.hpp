// Copyright (c) 2022-2024 Simons Foundation
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You may obtain a copy of the License at
//     https://www.gnu.org/licenses/gpl-3.0.txt
//
// Authors: Nikita Kavokine, Hao Lu, Olivier Parcollet, Nils Wentzell

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

    /// Perpendicular spin-spin correlation function :math:`\langle s_x(\tau) s_x(0) \rangle`.
    std::optional<gf<imtime>> sperp_tau;

    /// Density-density static correlation function :math:`\langle n_a(0) n_b(0) \rangle`.
    std::optional<nda::matrix<double>> nn_static;

    /// Density per color.
    nda::array<double, 1> densities;

    /// Delta perturbation order histogram
    std::optional<triqs::stat::histogram> pert_order_histo_Delta;

    /// J_perp perturbation order histogram
    std::optional<triqs::stat::histogram> pert_order_histo_Jperp;

    /// State histogram
    std::optional<nda::vector<double>> state_hist;

    /// Average sign
    double sign;
  };

  /// writes all containers to hdf5 file
  void h5_write(h5::group h5group, std::string subgroup_name, results_t const &c);

  /// reads all containers from hdf5 file
  void h5_read(h5::group h5group, std::string subgroup_name, results_t &c);

} // namespace triqs_ctseg
