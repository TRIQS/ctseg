#include "configuration.hpp"
#include "logs.hpp"
#include <iomanip>

// ===================  Functions to manipulate segments ===================

// Split a cyclic segment into 2 segment, attached to beta and 0
std::pair<segment_t, segment_t> split_cyclic_segment(segment_t const &s) {
  return {{tau_t::beta(), s.tau_cdag}, {s.tau_c, tau_t::zero()}};
}

// ---------------------------

// Check whether time is in segment [tau_c,tau_cdag[.
bool tau_in_seg(tau_t const &tau, segment_t const &seg) {
  // if the segment is cyclic, we check if tau in the 2 parts
  if (is_cyclic(seg)) {
    auto [sl, sr] = split_cyclic_segment(seg);
    return tau_in_seg(tau, sl) or tau_in_seg(tau, sr);
  }
  return (tau <= seg.tau_c and tau > seg.tau_cdag); // ! decreasing time order
}

// ---------------------------

// Overlap between two (possibly cyclic) segments.
double overlap(segment_t const &s1, segment_t const &s2) {
  // first treat the cyclic case
  if (is_cyclic(s1)) {
    auto [sl, sr] = split_cyclic_segment(s1);
    return overlap(sl, s2) + overlap(sr, s2);
    //return overlap(segment_t{tau_t::beta(), s1.tau_cdag}, s2) + overlap(segment_t{s1.tau_c, tau_t::zero()}, s2);
  }
  if (is_cyclic(s2)) return overlap(s2, s1); // s1 is not cyclic any more
  //return overlap(s1, segment_t{tau_t::beta(), s2.tau_cdag}) + overlap(s1, segment_t{s2.tau_c, tau_t::zero()});
  if (s1.tau_cdag >= s2.tau_c or s2.tau_cdag >= s1.tau_c) return 0;
  // last case
  tau_t tau_start = std::min(s1.tau_c, s2.tau_c);
  tau_t tau_end   = std::max(s1.tau_cdag, s2.tau_cdag);
  return double(tau_start - tau_end);
};

// =================== Functions to manipulate std::vector<segment_t> ========

// Find the closest segment on the left of seg.
// If there is none, returns the first on the right
// the list shoud not be empty
vec_seg_iter_t find_segment_left(std::vector<segment_t> const &seglist, segment_t const &seg) {
  auto seg_iter = std::upper_bound(seglist.begin(), seglist.end(), seg);
  return (seg_iter == seglist.begin()) ? seg_iter : --seg_iter;
}
// ---------------------------

//// Length occupied by all segments for a given color
//double density(std::vector<segment_t> const &seglist) {
//if (seglist.empty()) return 0;
//double result = 0;
//for (auto const &seg : seglist) result += double(seg.tau_c - seg.tau_cdag);
//return result;
//}
// ---------------------------

// Find density in seglist to the right of time tau.
int n_tau(tau_t const &tau, std::vector<segment_t> const &seglist) {
  if (seglist.empty()) return 0;
  auto it = find_segment_left(seglist, segment_t{tau, tau});
  return (tau_in_seg(tau, *it) or tau_in_seg(tau, seglist.back())) ? 1 : 0;
}

// ---------------------------
#if 0
// Flip seglist
std::vector<segment_t> flip(std::vector<segment_t> const &sl) {
  if (sl.empty()) // Flipped seglist is full line
    return {segment_t{tau_t::beta(), tau_t::zero()}};
  if (sl.size() == 1 and is_full_line(sl[0])) // Do nothing: flipped config empty
    return {};
  // else Swap c and cdag
  auto fsl = std::vector<segment_t>{};
  fsl.reserve(sl.size() + 1);
  for (auto i : range(sl.size())) {
    if (is_cyclic(sl.back())) {
      long ind = (i == 0) ? long(sl.size()) - 1 : i - 1;
      fsl.emplace_back(segment_t{sl[ind].tau_cdag, sl[i].tau_c, sl[ind].J_cdag, sl[i].J_c});
    } else {
      long ind = (i == sl.size() - 1) ? 0 : i + 1;
      fsl.emplace_back(segment_t{sl[i].tau_cdag, sl[ind].tau_c, sl[i].J_cdag, sl[ind].J_c});
    }
  } // loop over segs
  return fsl;
}
#endif

// Flip seglist
// FIXME : Recheck and merge
// do not break the if loop with the condition ...
std::vector<segment_t> flip(std::vector<segment_t> const &sl) {
  if (sl.empty()) // Flipped seglist is full line
    return {segment_t::full_line()};

  if (sl.size() == 1 and is_full_line(sl[0])) // Do nothing: flipped config empty
    return {};

  long N   = sl.size();
  auto fsl = std::vector<segment_t>(N); // NB must be () here, not {} !
  if (is_cyclic(sl.back()))
    for (auto i : range(N)) {
      long ind = (i == 0) ? N - 1 : i - 1;
      fsl[i]   = segment_t{sl[ind].tau_cdag, sl[i].tau_c, sl[ind].J_cdag, sl[i].J_c};
    }
  else
    for (auto i : range(N)) {
      long ind = (i == N - 1) ? 0 : i + 1;
      fsl[i]   = segment_t{sl[i].tau_cdag, sl[ind].tau_c, sl[i].J_cdag, sl[ind].J_c};
    }
  return fsl;
}

// ---------------------------

// Overlap between segment and a list of segments.
double overlap(std::vector<segment_t> const &seglist, segment_t const &seg) {
  if (seglist.empty()) return 0;
  // If seg is cyclic, need to split it because of the condition in the for later
  if (is_cyclic(seg)) {
    auto [sl, sr] = split_cyclic_segment(seg);
    return overlap(seglist, sl) + overlap(seglist, sr);
  }

  // Compute overlap with all segments
  double result = 0;
  // R is the iterator on the segment strictly after seg or end.
  // L the segment before or begin
  auto R    = std::upper_bound(seglist.begin(), seglist.end(), seg);
  auto L    = (R == seglist.begin()) ? R : R - 1;
  auto last = seglist.end() - 1;
  for (auto it = L; it != last and it->tau_c > seg.tau_cdag; ++it) //
    result += overlap(*it, seg);

  result += overlap(*last, seg); // the last can be cyclic, hence be unreached due to it->tau_c condition
  // nb : overlap is ok to call on cyclic segment
  return result;
}

// ---------------------------

// Checks if two segments are completely disjoint. s2 might be cyclic, not s1
inline bool disjoint(segment_t const &s1, segment_t const &s2) {
  if (is_cyclic(s2)) {
    auto [sl, sr] = split_cyclic_segment(s2);
    return disjoint(s1, sl) and disjoint(s1, sr);
  }
  return s1.tau_cdag > s2.tau_c or s2.tau_cdag > s1.tau_c; // symmetric s1 s2
}
// ---------------------------
#if 0
// Checks if segment is movable to a given color
bool is_insertable(std::vector<segment_t> const &seglist, segment_t const &seg) {
  if (seglist.empty()) return true;
  // If seg is cyclic, split it
  if (is_cyclic(seg)) {
    auto [sl, sr] = split_cyclic_segment(seg);
    return is_insertable(seglist, sl) and is_insertable(seglist, sr);
  }
  bool ok = true;
  // In case last segment in list is cyclic, split it and check its overlap with seg
  if (is_cyclic(seglist.back())) {
    auto [sl, sr] = split_cyclic_segment(seglist.back());
    ok            = disjoint(seg, sl) and disjoint(seg, sr);
  } else
    ok = disjoint(seg, seglist.back());
  if (not ok) return false;

  // Check overlap of seg with the remainder of seglist
  for (auto it = find_segment_left(seglist, seg); it->tau_c >= seg.tau_cdag and it != --seglist.end(); ++it) {
    if (not disjoint(seg, *it)) return false;
  }
  return true;
}
#endif

// ---------------------------

// Checks if segment is movable to a given color
bool is_insertable_into(segment_t const &seg, std::vector<segment_t> const &seglist) {
  if (seglist.empty()) return true;

  // If seg is cyclic, split it
  if (is_cyclic(seg)) {
    auto [sl, sr] = split_cyclic_segment(seg);
    return is_insertable_into(sl, seglist) and is_insertable_into(sr, seglist);
  }

  // R is the iterator on the segment strictly after seg or end.
  // L the segment before or begin
  // Then the segment is insertable iif it does not overlap with L not with R (if not end)
  // Proof : it overlaps with any segment before of equal L iff it does with L
  //         it overlaps with any segment after of equal L iif it does with R
  // see all cases.
  // 1-  L------       R-----  : s in [L,R]
  //           s------
  auto R = std::upper_bound(seglist.begin(), seglist.end(), seg);
  auto L = (R == seglist.begin()) ? R : R - 1;
  if (not disjoint(seg, *L)) return false;
  if (R != seglist.end() and not disjoint(seg, *R)) return false;
  // We must recheck the last segment as it may be cyclic (it might also have been R, in which case it is superfluous but ok)
  if (not disjoint(seg, seglist.back())) return false;
  return true;
}

// ---------------------------

// FIXME : why not pass a setgment ??
// FIXME : why not pass a view of K
// Contribution of the dynamical interaction kernel K to the overlap between a segment and a list of segments.
double K_overlap(std::vector<segment_t> const &seglist, tau_t const &tau_c, tau_t const &tau_cdag,
                 gf<imtime, matrix_valued> const &K, int c1, int c2) {

  // seglist empty covered by the loop
  double result = 0;
  for (auto const &s : seglist) {
    result += real(K(double(tau_c - s.tau_c))(c1, c2) + K(double(tau_cdag - s.tau_cdag))(c1, c2)
                   - K(double(tau_cdag - s.tau_c))(c1, c2) - K(double(tau_c - s.tau_cdag))(c1, c2));
  }
  return result;
}

// ---------------------------

// FIXME : the is_c is a HACK !
// Contribution of the dynamical interaction kernel K to the overlap between an operator and a list of segments.
double K_overlap(std::vector<segment_t> const &seglist, tau_t const &tau, bool is_c, gf<imtime, matrix_valued> const &K,
                 int c1, int c2) {
  double result = 0;
  // The order of the times is important for the measure of F
  for (auto const &s : seglist) {
    result += real(K(double(s.tau_c - tau))(c1, c2) - K(double(s.tau_cdag - tau))(c1, c2));
  }
  return is_c ? result : -result;
}

// ===================  Functions to manipulate config ===================

int n_at_boundary(configuration_t const &config, int color) {
  auto const &sl = config.seglists[color];
  if (sl.empty()) return 0;
  return (is_cyclic(sl.back()) or is_full_line(sl.back())) ? 1 : 0;
}

// ---------------------------

// Find segments corresponding to bosonic line
std::pair<vec_seg_iter_t, vec_seg_iter_t> find_spin_segments(int line_idx, configuration_t const &config) {
  auto const &line    = config.Jperp_list[line_idx];
  auto const &sl_up   = config.seglists[0];
  auto const &sl_down = config.seglists[1];
  // In spin up line, the c conneted to the J line is a at tau_Sminus
  auto c_up  = segment_t{line.tau_Sminus, line.tau_Sminus};
  auto it_up = std::lower_bound(sl_up.cbegin(), sl_up.cend(), c_up);
  // In spin down line, the c conneted to the J line is a at tau_Splus
  auto c_down  = segment_t{line.tau_Splus, line.tau_Splus};
  auto it_down = std::lower_bound(sl_down.cbegin(), sl_down.cend(), c_down);
  return {it_up, it_down};
}

// ---------------------------

// Sign of a configuration
double config_sign(configuration_t const &config, std::vector<det_t> const &dets) {
  double sign = 1.0;
  // For every color we compute the sign of the permutation that takes
  // [(c_dag c) (c_dag c) (c_dag c) ...] with the cdag in increasing time order
  // to the completely time-ordered list of operators (with increasing time)
  for (auto c : range(config.n_color())) {
    auto s = long(dets[c].size());
    if (s != 0) {
      // We first compute the sign of the permutation that takes
      // [(c_dag c) (c_dag c) (c_dag c) ...] with the cdag time-ordered to
      // [(c_dag c_dag ... cdag)(c c ... c)] with the c and c_dag time-ordered
      if ((s * (s - 1) / 2) % 2 == 1) sign *= -1;
      // We then compute the sign of the permutation that takes
      // [(c_dag c_dag ... cdag)(c c ... c)] with the c and c_dag time-ordered
      // to the completely time-ordered list of operators
      int idx_c = 0, idx_cdag = 0;
      for (int n = 0; n < 2 * s - 1; ++n) {
        if (dets[c].get_x(idx_cdag).first < dets[c].get_y(idx_c).first) {
          if (idx_cdag < s - 1) ++idx_cdag;
        } else {
          // Count the number of transpositions
          if ((s - idx_cdag) % 2 == 1) sign *= -1;
          if (idx_c < s - 1) ++idx_c;
        }
      }
    }
  }
  return sign;

  /*   
    // Old sign computation, works only without J_perp
    for (auto const &[c, sl] : itertools::enumerate(config.seglists)) {
    if (not sl.empty()) {
      bool starts_with_dagger = false;
      auto s                  = dets[c].size();
      if (s != 0) starts_with_dagger = dets[c].get_x(s - 1) > dets[c].get_y(s - 1);
      if (starts_with_dagger) sign *= (s % 2 == 0) ? 1 : -1; 

    }
  } */
}

// ---------------------------

// Find the indices of the segments whose cdag are in ]wtau_left,wtau_right[
std::vector<long> cdag_in_window(tau_t const &wtau_left, tau_t const &wtau_right,
                                 std::vector<segment_t> const &seglist) {
  if (seglist.empty()) return {}; // should never happen, but protect

  if (wtau_left < wtau_right) {
    auto left_list  = cdag_in_window(tau_t::beta(), wtau_right, seglist);
    auto right_list = cdag_in_window(wtau_left, tau_t::zero(), seglist);
    // concatenate
    left_list.insert(left_list.end(), right_list.begin(), right_list.end());
    return left_list;
  }
  // FIXME : REREAD ???
  std::vector<long> found_indices;
  found_indices.reserve(seglist.size());
  for (auto it = find_segment_left(seglist, segment_t{wtau_left, wtau_left});
       it->tau_cdag > wtau_right and it != --seglist.end(); ++it) {
    if (it->tau_cdag < wtau_left) found_indices.push_back(std::distance(seglist.cbegin(), it));
  }
  // Check separately for last segment (may be cyclic)
  if (seglist.back().tau_cdag < wtau_left and seglist.back().tau_cdag > wtau_right)
    found_indices.push_back(seglist.size() - 1);
  return found_indices;
}

// FIXME : do we have TESTS ???

/*
// Simpler implementation for test purposes

// Find the indices of the segments whose cdag are in ]wtau_left,wtau_right[
std::vector<long> cdag_in_window(tau_t const &wtau_left, tau_t const &wtau_right,
                                 std::vector<segment_t> const &seglist) {
  std::vector<long> found_indices;
  if (wtau_left < wtau_right) {
    auto left_list  = cdag_in_window(tau_t::beta(), wtau_right, seglist);
    auto right_list = cdag_in_window(wtau_left, tau_t::zero(), seglist);
    found_indices   = left_list;
    for (auto const &[i, idx] : itertools::enumerate(right_list)) found_indices.push_back(right_list[i]);
    return found_indices;
  }
  if (seglist.empty()) return found_indices; // should never happen, but protect
  found_indices.reserve(seglist.size());
  for (auto const &[i, seg] : itertools::enumerate(seglist)) {
    if (seg.tau_cdag < wtau_left and seg.tau_cdag > wtau_right) found_indices.push_back(i);
  }
  return found_indices;
}
*/

// ---------------------------

// Print config
std::ostream &operator<<(std::ostream &out, configuration_t const &config) {
  for (auto const &[c, sl] : itertools::enumerate(config.seglists)) {
    out << '\n';
    for (auto const &[i, seg] : itertools::enumerate(sl))
      out << "Color " << c << ". Position " << i << " : [ J:" << seg.J_c << " " << seg.tau_c << ", " << seg.tau_cdag
          << " J:" << seg.J_cdag << "]\n";
  }
  out << "\nSpin lines : \n";
  for (auto const &[i, line] : itertools::enumerate(config.Jperp_list)) {
    out << "S_minus : [" << line.tau_Sminus << "] S_plus : [" << line.tau_Splus << "]\n";
  }
  return out;
}
