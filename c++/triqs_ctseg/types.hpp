#pragma once
#include <fstream>
#include <triqs/gfs.hpp>
#include <triqs/hilbert_space/fundamental_operator_set.hpp>
#include <triqs/mc_tools/random_generator.hpp>
#include <triqs/operators/many_body_operator.hpp>
#include <triqs/operators/util/extractors.hpp>
#include <triqs/stat.hpp>
#include <triqs_ctseg/tau_t.hpp>
#include <triqs/det_manip.hpp>

#include "./tau_t.hpp"

using namespace triqs::gfs;
using namespace triqs::mesh;
using namespace triqs::stat;
using namespace triqs::utility;
using namespace triqs::operators;
using namespace triqs::hilbert_space;
using namespace h5;
using Op = many_body_operator;

// One-particle Green's function types
//using G_tau_t = block_gf<imtime, matrix_valued>;
//using G_iw_t  = block_gf<imfreq, matrix_valued>;

/// A lambda to adapt the Delta function for the call of the det.
/*struct delta_block_adaptor {*/
  //gf<imtime, matrix_real_valued> delta; // make a copy. Needed in the real case anyway.

  //double operator()(std::pair<tau_t, int> const &x, std::pair<tau_t, int> const &y) const {
    //double res = delta(double(x.first - y.first))(x.second, y.second);
    //return (x.first >= y.first ? res : -res); // x,y first are time_pt, wrapping is automatic in
                                              //// the - operation, but need to compute the sign
  //}
/*};*/

//using det_t = triqs::det_manip::det_manip<delta_block_adaptor>;
