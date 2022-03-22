#pragma once

#include "params.hpp"
#include <triqs/gfs.hpp>

// Group the inputs of the solvers.

struct inputs_t {

  block_gf<imtime> delta;
  gf<imtime, matrix_valued> jperpt; // perp spin retarded interaction kernel
  gf<imtime, matrix_valued> d0t;
};
