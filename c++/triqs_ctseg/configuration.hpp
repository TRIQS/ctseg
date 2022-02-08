#pragma once
#include <vector>

/// Segment: (time of c, time of c^dagger)
struct segment_t {
  qmc_time_t tau_c, tau_cdag; // time of c and cdag
  bool J_at_start = false, J_at_end = false;
};

// segment can be cyclic : tau_cdag > t_c : must be the last one ... FIXME : check invariant ...

// We order the segments by decreasing tau of the c operator
inline bool operator<(segment_t const &s1, segment_t const &s2) { return s1.tau_c > s2.tau_c; };


struct jperp_line_t {
  qmc_time_t tau_Splus, tau_Sminus; // times of the S+, S-
  //bool Splus_at_left;               /// USEFUL ???
};

struct configuration_t {
  std::vector<std::vector<segment_t>> seglists; // list of segment per color : seglist[color] is ORDERED on tau, with decreasing order.
  std::vector<jperp_line_t> Jperp_list;
};

// --- functions to manipulate config ---

inline bool operator<(segment_t const &s1, segment_t const &s2) { return s1.tau_c > s2.tau_c; }; // returns 1 if s1 is left of s2 

// returns index of first segment to the left of tau
inline auto find_segment_left(std::vector<segment_t> const &seglist, segment_t const &seg) {
  auto seg_index = std::upper_bound(seglist.begin(), seglist.end(), seg);
  return seg_index--; 
};

inline double overlap(std::vector<segment_t> const &seglist, segment_t const &seg) {
  auto ind = find_segment_left(seglist, seg);
  double result = 0;
  while (seglist[ind].tau_c < seg.tau_cdag) {
    

  }
  
  return result;
};


// --------- DEBUG code --------------
// print config + h5 config

inline void check_invariant(std::vector<segment_t> const &seglist) {
  // debug mode : check ordered.

  // position of J, etc...
}


