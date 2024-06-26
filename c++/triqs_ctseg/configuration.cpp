#include "configuration.hpp"

namespace triqs_ctseg {

  // ===================  Functions to manipulate segments ===================

  // Split a cyclic segment into 2 segment, attached to beta and 0
  std::pair<segment_t, segment_t> split_cyclic_segment(segment_t const &s) {
    return {{tau_t::beta(), s.tau_cdag}, {s.tau_c, tau_t::zero()}};
  }

  // -------------- internal -------------

  // Checks if two segments are completely disjoint.
  // !! s2 might be cyclic, not s1
  bool disjoint(segment_t const &s1, segment_t const &s2) {
    if (is_cyclic(s2)) {
      auto [sl, sr] = split_cyclic_segment(s2);
      return disjoint(s1, sl) and disjoint(s1, sr);
    }
    return s1.tau_cdag > s2.tau_c or s2.tau_cdag > s1.tau_c; // symmetric s1 s2
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
    }
    if (is_cyclic(s2)) return overlap(s2, s1); // s1 is not cyclic any more
    // [c1 cd1]                            [c1 cd1]
    //          [c2 cd2]   OR    [c2 cd2]
    if (s1.tau_cdag >= s2.tau_c or s2.tau_cdag >= s1.tau_c) return 0;
    // last case
    // [c1     cd1]                   [c1 cd1]
    //     [c2      cd2]      [c2            cd2]
    tau_t tau_start = std::min(s1.tau_c, s2.tau_c);       // ! most RIGHT as ordered by >  !
    tau_t tau_end   = std::max(s1.tau_cdag, s2.tau_cdag); // most LEFT
    assert(tau_start >= tau_end);
    return double(tau_start - tau_end);
  };

  // =================== Functions to manipulate std::vector<segment_t> ========

  vec_seg_iter_t lower_bound(std::vector<segment_t> const &seglist, tau_t const &tau) {
    // comparison is s.tau > tau as in the tau_t comparison
    return std::lower_bound(seglist.begin(), seglist.end(), tau, [](auto &&s, auto &&t) { return s.tau_c > t; });
  }

  // ------------ internal ---------------

  // Iterator on the closest segment on the left of seg.
  // If there is none, returns the first on the right (or end)
  // the list shoud not be empty
  vec_seg_iter_t find_segment_left(std::vector<segment_t> const &seglist, segment_t const &seg) {
    auto seg_iter = std::upper_bound(seglist.begin(), seglist.end(), seg);
    return (seg_iter == seglist.begin()) ? seg_iter : --seg_iter;
  }

  // ---------------------------

  int n_at_boundary(std::vector<segment_t> const &sl) {
    if (sl.empty()) return 0;
    return (is_cyclic(sl.back()) or is_full_line(sl.back())) ? 1 : 0;
  }

  // ---------------------------

  // Find density in seglist to the right of time tau.
  int n_tau(tau_t const &tau, std::vector<segment_t> const &seglist) {
    if (seglist.empty()) return 0;
    auto it = find_segment_left(seglist, segment_t{tau, tau});
    return (tau_in_seg(tau, *it) or tau_in_seg(tau, seglist.back())) ? 1 : 0;
  }

  // ---------------------------
  // Flip seglist
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
    // first loop on all segment but the last one
    auto last = seglist.end() - 1;
    for (auto it = find_segment_left(seglist, seg); it != last and it->tau_c > seg.tau_cdag; ++it) //
      result += overlap(*it, seg);

    // the last can be cyclic, hence be unreached due to it->tau_c condition
    // nb : overlap is ok to call on cyclic segment
    result += overlap(*last, seg);
    return result;
  }

  // ---------------------------

  // Checks if segment is insertable to a given color
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
  // FIXME : do we have TESTS ???
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
    std::vector<long> found_indices;
    found_indices.reserve(seglist.size());
    auto last = seglist.end() - 1;
    auto it   = find_segment_left(seglist, segment_t{wtau_left, wtau_left});
    for (; it->tau_cdag > wtau_right and it != last; ++it)
      if (it->tau_cdag < wtau_left) found_indices.push_back(std::distance(seglist.cbegin(), it));

    // Check separately for last segment (may be cyclic)
    if (seglist.back().tau_cdag < wtau_left and seglist.back().tau_cdag > wtau_right)
      found_indices.push_back(seglist.size() - 1);
    return found_indices;
  }

  // ---------------------------

  // Contribution of the dynamical interaction kernel K to the overlap between a segment and a list of segments.
  // Computes the sum of the s_a s_b K(tau_a - tau_b) where s_a is 1 for cdag and - 1 for c
  double K_overlap(std::vector<segment_t> const &seglist, tau_t const &tau_c, tau_t const &tau_cdag,
                   gf<imtime, matrix_valued> const &K, int c1, int c2) {

    auto Ks = slice_target_to_scalar(K, c1, c2);

    // seglist empty covered by the loop
    double result = 0;
    for (auto const &s : seglist) {
      result += real(Ks(double(tau_c - s.tau_c)) + Ks(double(tau_cdag - s.tau_cdag)) - Ks(double(tau_cdag - s.tau_c))
                     - Ks(double(tau_c - s.tau_cdag)));
    }
    return result;
  }

  // ---------------------------

  // Contribution of the dynamical interaction kernel K to the overlap between an operator and a list of segments.
  // Computes the sum of the s_a s_b K(tau_a - tau_b) where s_a is 1 for cdag and - 1 for c
  double K_overlap(std::vector<segment_t> const &seglist, tau_t const &tau, bool is_c,
                   gf<imtime, matrix_valued> const &K, int c1, int c2) {
    auto Ks = slice_target_to_scalar(K, c1, c2);

    double result = 0;
    // The order of the times is important for the measure of F
    for (auto const &s : seglist) { result += real(Ks(double(s.tau_c - tau)) - Ks(double(s.tau_cdag - tau))); }
    return is_c ? result : -result;
  }

  // ---------------------------

  // List of operators containing all colors.
  // Time are ordered in decreasing order, in agreement with the whole physic literature.
  std::vector<colored_ops_t> colored_ordered_ops(std::vector<std::vector<segment_t>> const &seglists) {
    int c = 0;                           // index of color
    std::vector<colored_ops_t> ops_list; // list of all the operators
    for (auto const &seglist : seglists) {
      for (auto const &s : seglist) {
        ops_list.push_back({s.tau_c, c, false});
        ops_list.push_back({s.tau_cdag, c, true});
        if (is_cyclic(s)) {
          ops_list.push_back({tau_t::beta(), c, false});
          ops_list.push_back({tau_t::zero(), c, true});
        }
      }
      ++c;
    }
    std::sort(ops_list.begin(), ops_list.end(), [](const colored_ops_t &a, const colored_ops_t &b) {
      return b.tau < a.tau; // Note the order of b and a for descending sort
    });
    return ops_list;
  }

  // ===================  PRINTING ========================

  std::ostream &operator<<(std::ostream &out, std::vector<segment_t> const &sl) {
    out << '\n';
    for (auto const &[i, seg] : itertools::enumerate(sl))
      out << ". Position " << i << " : [ J:" << seg.J_c << " " << seg.tau_c << ", " << seg.tau_cdag
          << " J:" << seg.J_cdag << "]\n";
    return out;
  }

  // ---------------------------

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

  // ---------------------------

  std::ostream &operator<<(std::ostream &out, std::vector<colored_ops_t> const &col) {
    for (auto const &[i, co] : itertools::enumerate(col))
      out << "\n"
          << "i: " << i << ", tau: " << co.tau << ", color: " << co.color << ", " << (co.is_cdag ? "cdag" : "c");
    return out;
  }

} // namespace triqs_ctseg
