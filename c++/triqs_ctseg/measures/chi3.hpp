// Copyright (c) 2024--present, The Simons Foundation
// This file is part of TRIQS/ctseg and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#pragma once
#include "../configuration.hpp"
#include "../work_data.hpp"
#include "../results.hpp"

#ifdef TRIQS_WITH_NFFT
#include <triqs/utility/nfft/buffer.hpp>
#endif

namespace triqs_ctseg::measures {

#ifdef TRIQS_WITH_NFFT

  struct chi3 {

    chi3(params_t const &params, work_data_t const &wdata, configuration_t const &config, results_t &results);

    // chi3 is uncopyable due to nfft::buffer_t
    chi3(chi3 const &)            = delete;
    chi3(chi3 &&)                 = default;
    ~chi3()                       = default;
    chi3 &operator=(chi3 const &) = delete;
    chi3 &operator=(chi3 &&)      = delete;

    void accumulate(double s);
    void collect_results(mpi::communicator const &c);

    private:
    work_data_t const &wdata;
    configuration_t const &config;
    results_t &results;
    double beta;

    // DLR2D mesh for chi3 output
    triqs::mesh::dlr2d_imfreq chi3_mesh;

    // NFFT buffers for M: buf_arr(bl)(a, b) — one per orbital pair per block
    nda::array<nda::array<triqs::utility::nfft::buffer_t<2>, 2>, 1> buf_arr;

    // Intermediate M on DLR2D mesh (block_gf, matrix_valued)
    block_gf<dlr2d_imfreq, matrix_valued> M;

    // Target Matsubara frequencies for type3 NFFT
    std::vector<std::array<triqs::mesh::matsubara_freq, 2>> target_mf;

    // Per-color nw on bosonic mesh
    std::vector<gf<imfreq, scalar_valued>> nw;
    triqs::mesh::imfreq nw_mesh;

    // Output accumulator (vector of vectors, later packed into block2_gf)
    std::vector<std::vector<gf<dlr2d_imfreq, tensor_valued<4>>>> chi3_acc;

    // Block names for final output
    std::vector<std::string> block_names;

    double Z = 0;
  };

#endif // TRIQS_WITH_NFFT

} // namespace triqs_ctseg::measures
