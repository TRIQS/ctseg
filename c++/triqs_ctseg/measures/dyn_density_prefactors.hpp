// Copyright (c) 2026--present, The Simons Foundation
// Copyright (c) 2026--present, Max Planck Institute for Polymer Research, Mainz, Germany
// This file is part of TRIQS/ctseg and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#pragma once

#include "../configuration.hpp"
#include "../work_data.hpp"

namespace triqs_ctseg::measures {

  inline double static_prefactor(work_data_t const &wdata, configuration_t const &config, int color, tau_t tau) {
    double value = 0;
    for (auto const &[c, sl] : itertools::enumerate(config.seglists)) {
      if (c != color) value += wdata.U(c, color) * n_tau(tau, sl);
    }
    return value;
  }

  inline double retarded_density_prefactor(work_data_t const &wdata, configuration_t const &config, int color,
                                           tau_t tau) {
    if (not wdata.has_Dt) return 0;

    double value = 0;
    for (auto const &[c, sl] : itertools::enumerate(config.seglists))
      value -= K_overlap(sl, tau, false, wdata.Kprime, c, color);
    return value;
  }

  inline double retarded_density_jump(work_data_t const &wdata, int color) {
    return wdata.has_Dt ? -2 * real(wdata.Kprime(0)(color, color)) : 0;
  }

  inline double retarded_density_prefactor_integral(work_data_t const &wdata, configuration_t const &config, int color,
                                                    tau_t tau_left, tau_t tau_right) {
    if (not wdata.has_Dt) return 0;

    double value = 0;
    for (auto const &[c, sl] : itertools::enumerate(config.seglists)) {
      auto Ks = slice_target_to_scalar(wdata.K, c, color);
      for (auto const &seg : sl) {
        value += real(Ks(double(seg.tau_c - tau_right)) - Ks(double(seg.tau_c - tau_left))
                      - Ks(double(seg.tau_cdag - tau_right)) + Ks(double(seg.tau_cdag - tau_left)));
      }
    }
    return value;
  }

} // namespace triqs_ctseg::measures
