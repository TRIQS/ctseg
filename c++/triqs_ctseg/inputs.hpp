#pragma once

#include "params.hpp"
#include <triqs/gfs.hpp>

// Group the inputs of the solver.

struct inputs_t {

  block_gf<imtime> delta; // hybridization function 
  gf<imtime, matrix_valued> jperpt; // perpendicular spin-spin interaction
  block2_gf<imtime> d0t; // retarded density-density interaction
  
};

// h5_read/write
void h5_write(h5::group h5group, std::string subgroup_name, inputs_t const &s);
void h5_read(h5::group h5group, std::string subgroup_name, inputs_t &s);
