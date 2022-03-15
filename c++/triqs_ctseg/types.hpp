#pragma once
#include <fstream>
#include <triqs/gfs.hpp>
#include <triqs/hilbert_space/fundamental_operator_set.hpp>
#include <triqs/mc_tools/random_generator.hpp>
//#include <triqs/mesh.hpp>
#include <triqs/operators/many_body_operator.hpp>
#include <triqs/operators/util/extractors.hpp>
#include <triqs/stat.hpp>
#include <triqs/utility/time_pt.hpp>

// FIXME
#define EXT_DEBUG

using namespace triqs::gfs;
using namespace triqs::mesh;
using namespace triqs::stat;
using namespace triqs::utility;
using namespace triqs::operators;
using namespace triqs::hilbert_space;
using namespace h5;
using Op = many_body_operator;
using std::cos;
using std::cout;
using std::endl;
using std::sin;

// One-particle Green's function types
using G_tau_t            = block_gf<imtime, matrix_valued>;
using G_iw_t             = block_gf<imfreq, matrix_valued>;
using qmc_time_t         = triqs::utility::time_pt;
using qmc_time_factory_t = triqs::utility::time_segment;
