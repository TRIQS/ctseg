#include "configuration.hpp"

// ------------------- Invariants ---------------------------

void check_invariant(configuration_t const &config) {
  for (auto const &[c, sl] : itertools::enumerate(config.seglists))
    for (int i = 0; i < sl.size() - 1; ++i) {
      ALWAYS_EXPECTS((sl[i].tau_cdag > sl[i + 1].tau_c),
                     "Time order error between segment at position {} in config \n{}", i, config);
      ALWAYS_EXPECTS(not is_cyclic(sl[i]), "Segment at position {} should not by cyclic in config \n{}", i,
                     config); // only last segment can be cyclic
    }
}
// ---------------------------

// Make a list of time ordered (decreasing) operators

std::vector<std::tuple<qmc_time_t, int, bool>> make_time_ordered_op_list(configuration_t const &config) {

  std::vector<std::tuple<qmc_time_t, int, bool>> result;

  result.reserve(config.seglists.size()
                 * config.seglists[0].size()); // optional : simple heuristics to reserve some space

  for (auto const &[color, seglist] : itertools::enumerate(config.seglists))
    for (auto const &seg : seglist) {
      result.emplace_back(seg.tau_c, color, false);
      result.emplace_back(seg.tau_cdag, color, true);
    }

  std::sort(std::begin(result), std::end(result), std::greater{});
  return result;
}

// ---------------------------

// Overlap between two non-cyclic segments.
double overlap_seg(segment_t const &seg1, segment_t const &seg2) {
  if (seg1.tau_cdag >= seg2.tau_c || seg2.tau_cdag >= seg1.tau_c)
    return 0;
  else {
    qmc_time_t tau_start = std::min(seg1.tau_c, seg2.tau_c);
    qmc_time_t tau_end   = std::max(seg1.tau_cdag, seg2.tau_cdag);
    return double(tau_start - tau_end);
  };
};
// ---------------------------

// Overlap between segment and a list of segments.
double overlap(std::vector<segment_t> const &seglist, segment_t const &seg, qmc_time_factory_t const &fac) {
  if (seglist.empty()) return 0;

  auto zero = fac.get_lower_pt();
  auto beta = fac.get_upper_pt();

  // If seg is cyclic, split it
  if (seg.tau_c < seg.tau_cdag)
    return overlap(seglist, segment_t{beta, seg.tau_cdag}, fac) + overlap(seglist, segment_t{seg.tau_c, zero}, fac);
  double result = 0;
  // Isolate last segment
  auto last_seg = seglist.back();
  // In case last segment is cyclic, split it and compute its overlap with seg
  if (last_seg.tau_c < last_seg.tau_cdag)
    result += overlap_seg(seg, segment_t{beta, last_seg.tau_cdag}) + overlap_seg(seg, segment_t{last_seg.tau_c, zero});
  else
    result += overlap_seg(seg, last_seg);

  // Compute overlap of seg with the remainder of seglist
  for (auto it = find_segment_left(seglist, seg); it->tau_c < seg.tau_cdag && it != --seglist.end(); ++it) //
    result += overlap_seg(*it, seg);
  return result;
};
// ---------------------------

// Contribution of the dynamical interaction kernel K to the overlap between a segment and a list of segments.
double K_overlap(std::vector<segment_t> const &seglist, segment_t const &seg, gf<imtime, scalar_valued> const &K) {
  if (seglist.empty()) return 0;
  double result = 0;
  for (auto seg_in_list : seglist) {
    result += real(K(double(seg.tau_c - seg_in_list.tau_c)) + K(double(seg.tau_cdag - seg_in_list.tau_cdag))
                   - K(double(seg.tau_cdag - seg_in_list.tau_c)) - K(double(seg.tau_c - seg_in_list.tau_cdag)));
  }
  return result;
}
// ---------------------------

// Length occupied by all segments for a given color
double density(std::vector<segment_t> const &seglist) {
  double result = 0;
  for (auto const &seg : seglist) result += double(seg.tau_c - seg.tau_cdag);
  return result;
}
// ---------------------------

/// Find the state at tau = 0 or beta
std::vector<bool> boundary_state(configuration_t const &config) {
  int N = config.n_color();
  std::vector<bool> res(N);
  for (auto const &[c, sl] : itertools::enumerate(config.seglists))
    res[c] = (sl.empty() ? false : is_cyclic(sl.back())); // FIXME  BUG full line
  return res;
}

// ---------------------------

std::ostream &operator<<(std::ostream &out, configuration_t const &config) {
  for (auto const &sl : config.seglists) {
    out << '\n';
    for (auto const &seg : sl) out << "[" << double(seg.tau_c) << ", " << double(seg.tau_cdag) << "]    ";
  }
  return out;
}
