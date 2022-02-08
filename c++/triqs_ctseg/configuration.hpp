#pragma once
#include <vector>
#include "params.hpp"

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

// ------------------- Functions to manipulate config --------------------------

// Comparison of segments. Returns 1 if s1 is left of s2. 
inline bool operator<(segment_t const &s1, segment_t const &s2) { return s1.tau_c > s2.tau_c; }; 

// Overlap between two non-cyclic segments. 
inline double overlap_seg(segment_t const &seg1, segment_t const &seg2) {
  if (seg1.tau_cdag >= seg2.tau_c || seg2.tau_cdag >= seg1.tau_c) return 0; 
  else {
    qmc_time_t tau_start = std::min(seg1.tau_c,seg2.tau_c);
    qmc_time_t tau_end = std::max(seg1.tau_cdag,seg2.tau_cdag);
    return double(tau_start-tau_end);
  };
};

// Find index of first segment to the left of seg.
inline auto find_segment_left(std::vector<segment_t> const &seglist, segment_t const &seg) {
  auto seg_index = std::upper_bound(seglist.begin(), seglist.end(), seg);
  return seg_index--; 
};

// Overlap between segment and a list of segments. 
inline double overlap(std::vector<segment_t> const &seglist, segment_t const &seg) {
  if seglist.empty() return 0;
  // If seg is cyclic, split it
  if (seg.tau_c < seg.tau_cdag) return overlap(seglist,segment_t{params.beta,tau_cdag}) + overlap(seglist,segment_t{tau_c,0});
  double result = 0;
  // Isolate last segment 
  auto last_seg = seglist.back();
  // In case last segment is cyclic, split it and compute its overlap with seg
  if (last_seg.tau_c < last_seg.tau_cdag) {
    result += overlap_seg(seg,segment_t{params.beta,last_seg.tau_cdag}) + overlap_seg(seg,segment_t{last_seg.tau_c,0});
  }
  else result+= overlap_seg(seg,last_seg);
  // Compute overlap of seg with the remainder of seglist
  auto ind = find_segment_left(seglist, seg);
  while (seglist[ind].tau_c < seg.tau_cdag && ind != seglist.end()) {
    result += overlap_seg(seglist[ind],seg);
    ++ind;
  };
  return result; 
};

// Length occupied by all segments for a given color
inline double density(std::vector<segment_t> const &seglist){
  double result = 0;
  for (int i : seglist){
    result += seglist[i].tau_c - seglist[i].tau_cdag;
  };
  return result; 
};
  
  



// --------- DEBUG code --------------
// print config + h5 config

inline void check_invariant(std::vector<segment_t> const &seglist) {
  // debug mode : check ordered.

  // position of J, etc...
};


