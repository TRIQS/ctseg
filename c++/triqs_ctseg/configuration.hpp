#pragma once
#include <vector>
#include "params.hpp"

// Segment: (time of c, time of c^dagger)
struct segment_t {
  qmc_time_t tau_c, tau_cdag; // time of c and cdag
  bool J_at_start = false, J_at_end = false;
};

// segment can be cyclic : tau_cdag > t_c : must be the last one ... FIXME : check invariant ...

struct jperp_line_t {
  qmc_time_t tau_Splus, tau_Sminus; // times of the S+, S-
  //bool Splus_at_left;               /// USEFUL ???
};

struct configuration_t {
  std::vector<std::vector<segment_t>> seglists; // list of segment per color : seglist[color] is ORDERED on tau, with decreasing order.
  std::vector<jperp_line_t> Jperp_list;
};

// ------------------- Functions to manipulate config --------------------------

// Comparison of segments. Returns 1 if s1 is left of s2 (we order segments by decreasing time). 
inline bool operator<(segment_t const &s1, segment_t const &s2) { return s1.tau_c > s2.tau_c; }; 

// Overlap between two non-cyclic segments. 
double overlap_seg(segment_t const &seg1, segment_t const &seg2);

// Find index of first segment to the left of seg.
auto find_segment_left(std::vector<segment_t> const &seglist, segment_t const &seg);

// Overlap between segment and a list of segments. 
double overlap(std::vector<segment_t> const &seglist, segment_t const &seg);

// Length occupied by all segments for a given color
double density(std::vector<segment_t> const &seglist);
  
// --------- DEBUG code --------------
// print config + h5 config

inline void check_invariant(std::vector<segment_t> const &seglist) {
  // debug mode : check ordered.

  // position of J, etc...
};


