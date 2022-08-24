#pragma once
#include <triqs/gfs.hpp>
#include <triqs/det_manip.hpp>

#include "./tau_t.hpp"

using namespace triqs::gfs;
using namespace triqs::mesh;

using delta_target_t = matrix_valued;
using delta_scalar_t = typename delta_target_t::scalar_t;

// -----------     dets    ----------------------------------------

/// A lambda to adapt the Delta function for the call of the det.
struct delta_block_adaptor {
  gf<imtime, delta_target_t> delta; // make a copy. Needed in the real case anyway.

  delta_scalar_t operator()(std::pair<tau_t, int> const &x, std::pair<tau_t, int> const &y) const {
    delta_scalar_t res = delta(double(x.first - y.first))(x.second, y.second);
    return (x.first >= y.first ? res : -res); // x,y first are time_pt, wrapping is automatic in
                                              // the - operation, but need to compute the sign
  }
};

using det_t = triqs::det_manip::det_manip<delta_block_adaptor>;
