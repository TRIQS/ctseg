#pragma once

#include "params.hpp"
#include <triqs/gfs.hpp>

// Group the inputs of the solvers.

struct inputs_t { 

  block_gf<imtime> delta;

  gf<imtime, matrix_valued> jperpt; // perp spin retarded interaction kernel
  gf<imtime, matrix_valued> d0t;
// blocks structure: (N_blocks)^2 (e.g if G is (up,dn), then (up|up, up|dn,
  // dn|up, dn|dn) )
  //block_gf<imtime> kt;  // retarded interaction kernel
  //block_gf<imtime> kprimet;

  //gf<imtime> kperpprimet; // derivative of retarded interaction kernel
  //gf<imfreq> jperpw;      // perp spin retarded interaction kernel
  //gf<imtime> chipmt;      // chi-plusminus spinspin susceptibility

};

