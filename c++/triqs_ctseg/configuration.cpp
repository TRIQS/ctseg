#include "configuration.hpp"

// Overlap between two non-cyclic segments. 
double overlap_seg(segment_t const &seg1, segment_t const &seg2) {
  if (seg1.tau_cdag >= seg2.tau_c || seg2.tau_cdag >= seg1.tau_c) return 0; 
  else {
    qmc_time_t tau_start = std::min(seg1.tau_c,seg2.tau_c);
    qmc_time_t tau_end = std::max(seg1.tau_cdag,seg2.tau_cdag);
    return double(tau_start-tau_end);
  };
};

// Find index of first segment to the left of seg.
auto find_segment_left(std::vector<segment_t> const &seglist, segment_t const &seg) {
  auto seg_index = std::upper_bound(seglist.begin(), seglist.end(), seg);
  return seg_index--; 
};

// Overlap between segment and a list of segments. 
double overlap(std::vector<segment_t> const &seglist, segment_t const &seg) {
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
double density(std::vector<segment_t> const &seglist){
  double result = 0;
  for (int i : seglist){
    result += seglist[i].tau_c - seglist[i].tau_cdag;
  };
  return result; 
};