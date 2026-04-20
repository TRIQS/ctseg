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

  struct chi4 {

    chi4(params_t const &params, work_data_t const &wdata, configuration_t const &config, results_t &results);

    // chi4 is uncopyable due to nfft::buffer_t
    chi4(chi4 const &)            = delete;
    chi4(chi4 &&)                 = default;
    ~chi4()                       = default;
    chi4 &operator=(chi4 const &) = delete;
    chi4 &operator=(chi4 &&)      = delete;

    void accumulate(double s);
    void collect_results(mpi::communicator const &c);

    private:
    work_data_t const &wdata;
    results_t &results;
    double beta;

    // Meshes
    triqs::mesh::imfreq mesh_bosonic;
    triqs::mesh::imfreq mesh_fermionic;
    triqs::mesh::imfreq aux_mesh; // enlarged fermionic mesh for M

    // NFFT buffers for M: buf_arr(bl)(a, b) — one per orbital pair per block
    nda::array<nda::array<triqs::utility::nfft::buffer_t<2>, 2>, 1> buf_arr;

    // Intermediate M on enlarged uniform mesh (block_gf, matrix_valued)
    block_gf<prod<imfreq, imfreq>, matrix_valued> M;

    // Output accumulator (vector of vectors, later packed into block2_gf)
    std::vector<std::vector<gf<prod<imfreq, imfreq, imfreq>, tensor_valued<4>>>> chi4_acc;

    // Block names for final output
    std::vector<std::string> block_names;

    double Z = 0;
  };

#endif // TRIQS_WITH_NFFT

} // namespace triqs_ctseg::measures
