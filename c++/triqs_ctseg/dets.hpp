#pragma once
#include <triqs/gfs.hpp>
#include <triqs/det_manip.hpp>

#include "./tau_t.hpp"

using namespace triqs::gfs;
using namespace triqs::mesh;

// -----------     dets    ----------------------------------------

/// A lambda to adapt the Delta function for the call of the det.
struct delta_block_adaptor {
  gf<imtime, matrix_real_valued> delta; // make a copy. Needed in the real case anyway.

  double operator()(std::pair<tau_t, int> const &x, std::pair<tau_t, int> const &y) const {
    double res = delta(double(x.first - y.first))(x.second, y.second);
    return (x.first >= y.first ? res : -res); // x,y first are time_pt, wrapping is automatic in
                                              // the - operation, but need to compute the sign
  }
};

using det_t = triqs::det_manip::det_manip<delta_block_adaptor>;


