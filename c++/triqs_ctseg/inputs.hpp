#pragma once

#include "params.hpp"
#include <triqs/gfs.hpp>

// Group the inputs of the solvers.

struct inputs_t {

  block_gf<imtime> delta;
  gf<imtime, matrix_valued> jperpt; // perp spin retarded interaction kernel
  gf<imtime, matrix_valued> d0t;
};

// h5_read/write
void h5_write(h5::group h5group, std::string subgroup_name, inputs_t const &s);
void h5_read(h5::group h5group, std::string subgroup_name, inputs_t &s);
