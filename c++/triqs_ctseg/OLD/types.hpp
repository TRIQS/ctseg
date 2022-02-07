/*******************************************************************************
 * CTSEG: TRIQS hybridization-expansion segment solver
 *
 * Copyright (C) 2013-2018 by T. Ayral, H. Hafermann, P. Delange, M. Ferrero, O.
 *Parcollet
 *
 * CTSEG is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * CTSEG is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * CTSEG. If not, see <http://www.gnu.org/licenses/>.
 ******************************************************************************/
#pragma once
#include <fstream>
#include <triqs/gfs.hpp>
#include <triqs/hilbert_space/fundamental_operator_set.hpp>
#include <triqs/mc_tools/random_generator.hpp>
#include <triqs/mesh.hpp>
#include <triqs/operators/many_body_operator.hpp>
#include <triqs/operators/util/extractors.hpp>
#include <triqs/stat.hpp>

namespace triqs_ctseg {

using namespace triqs::gfs;
using namespace triqs::mesh;
using namespace triqs::stat;
using namespace nda;
using namespace triqs::utility;
using namespace h5;
using std::cos;
using std::cout;
using std::endl;
using std::sin;

// Legendre Green function type
using g_l_t = block_gf<triqs::gfs::legendre, matrix_valued>;
using g_l_vt = block_gf_view<triqs::gfs::legendre, matrix_valued>;

// two-frequency container class
using block_f_om_nu_tv_t = block_gf<prod<imfreq, imfreq>, tensor_valued<3>>;
using block_f_om_nu_tv_vt =
    block_gf_view<prod<imfreq, imfreq>, tensor_valued<3>>;

// four-leg, three-frequency container class
using gf_3w_container_t =
    block_gf<prod<imfreq, imfreq, imfreq>, tensor_valued<4>>;

using namespace triqs::utility;
using namespace triqs::operators;
using namespace triqs::hilbert_space;
using Op = many_body_operator;
using gf_struct_t = triqs::hilbert_space::gf_struct_t;

} // namespace triqs_ctseg
