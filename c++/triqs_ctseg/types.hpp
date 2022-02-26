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

// spdlog
#ifdef EXT_DEBUG
#define SPDLOG_ACTIVE_LEVEL SPDLOG_LEVEL_TRACE
#else
#define SPDLOG_ACTIVE_LEVEL SPDLOG_LEVEL_OFF
#endif
#include "spdlog/spdlog.h"
#include "spdlog/fmt/ostr.h" 

#define LOG(...) SPDLOG_TRACE(__VA_ARGS__)

  // checked always, even in production
#define ALWAYS_EXPECTS(Condition, ErrorMessage, ...)                                                                                                        \
  if (not (Condition)) {                                                                                                                               \
    SPDLOG_CRITICAL(ErrorMessage, __VA_ARGS__);                                                                                                 \
    throw std::runtime_error("Assertion Error, cf log");                                                                                                 \
  }


// FIXME 

using namespace triqs::gfs;
using namespace triqs::mesh;
using namespace triqs::stat;
using namespace triqs::utility;
using namespace h5;
using std::cos;
using std::cout;
using std::endl;
using std::sin;

using namespace triqs::utility;
using namespace triqs::operators;
using namespace triqs::hilbert_space;
using Op = many_body_operator;


