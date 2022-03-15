#include "move_segment.hpp"
#include "../logs.hpp"

namespace moves {

  // FIXME WHY HERE ??
  // return positivie : do_overlap ?
  // Checks if two segments overlap (even just at their boundaries)
  bool move_segment::no_overlap(segment_t seg1, segment_t seg2) {
    if (seg1.tau_c >= seg2.tau_cdag && seg1.tau_c <= seg2.tau_c) return false;
    if (seg1.tau_cdag >= seg2.tau_cdag && seg1.tau_c <= seg2.tau_c) return false;
    return true;
  }

  // -------------------------------

  // Checks if a segment is movable to a color
  bool move_segment::is_movable(std::vector<segment_t> const &seglist, segment_t const &seg) {
    bool result = true;
    if (seglist.empty()) return result;
    // If seg is cyclic, split it
    if (seg.tau_c < seg.tau_cdag)
      return is_movable(seglist, segment_t{wdata.qmc_beta, seg.tau_cdag})
         && is_movable(seglist, segment_t{seg.tau_c, wdata.qmc_zero});
    // Isolate last segment
    segment_t last_seg = seglist.back();
    // In case last segment is cyclic, split it and check its overlap with seg
    if (last_seg.tau_c < last_seg.tau_cdag) {
      result = result && no_overlap(seg, segment_t{wdata.qmc_beta, last_seg.tau_cdag})
         && no_overlap(seg, segment_t{last_seg.tau_c, wdata.qmc_zero});
    } else
      result = result && no_overlap(seg, last_seg);
    // Check overlap of seg with the remainder of seglist
    for (auto it = find_segment_left(seglist, seg); it->tau_c < seg.tau_cdag && it != --seglist.end(); ++it) {
      result = result && no_overlap(*it, seg);
    }
    return result;
  }

  // -------------------------------

  double move_segment::attempt() {

    SPDLOG_TRACE("\n =================== ATTEMPT MOVE ================ \n");

    // ------------ Choice of segment and colors --------------
    // Select origin color
    origin_color = rng(wdata.n_color);
    auto &sl     = config.seglists[origin_color];
    SPDLOG_TRACE("Moving from color {}", origin_color);

    // If color has no segments, nothing to move
    if (sl.empty()) return 0;

    // Select segment to move
    origin_index                 = rng(sl.size());
    origin_segment               = sl[origin_index];
    auto qmc_length              = origin_segment.tau_c - origin_segment.tau_cdag;
    bool moving_full_line        = (qmc_length == wdata.qmc_beta);
    auto proposed_segment_length = double(qmc_length);

    SPDLOG_TRACE("Moving c at {}, cdag at {}", origin_segment.tau_c, origin_segment.tau_cdag);

    // Select destination color
    destination_color = rng(wdata.n_color - 1);
    if (destination_color >= origin_color) ++destination_color;
    auto &dsl = config.seglists[destination_color];

    // Reject if chosen segment overlaps with destination color
    if (not is_movable(dsl, origin_segment)) return 0;
    destination_it = ++find_segment_left(dsl, origin_segment);

    // ------------  Trace ratio  -------------
    double ln_trace_ratio = (wdata.mu(destination_color) - wdata.mu(origin_color)) * proposed_segment_length;
    double trace_ratio    = std::exp(ln_trace_ratio);

    // ------------  Det ratio  ---------------

    // pos is the position of the proposed segment if inserted, converted from iterator to int
    long dest_index = std::distance(destination_it, dsl.cbegin());
    // We insert tau_cdag as a line (first index) and tau_c as a column (second index). The index always corresponds to the
    // segment the tau_c/tau_cdag belongs to.
    double det_ratio = 1.0;
    if (not moving_full_line) {
      det_ratio = wdata.dets[destination_color].try_insert(dest_index, dest_index, {origin_segment.tau_cdag, 0},
                                                           {origin_segment.tau_c, 0})
         * wdata.dets[origin_color].try_remove(origin_index, origin_index);
    }

    // ------------  Proposition ratio ------------
    double prop_ratio = (int(dsl.size()) + 1) / int(sl.size());

    SPDLOG_TRACE("trace_ratio  = {}, prop_ratio = {}, det_ratio = {}", trace_ratio, prop_ratio, det_ratio);

    return trace_ratio * det_ratio * prop_ratio;
  }

  //--------------------------------------------------

  double move_segment::accept() {

    SPDLOG_TRACE("\n - - - - - ====> ACCEPT - - - - - - - - - - -\n");

    // Update the dets
    wdata.dets[origin_color].complete_operation();
    wdata.dets[destination_color].complete_operation();

    double sign_ratio = 1;

    // Add the segment at destination and update sign
    auto &dsl = config.seglists[destination_color];
    if (is_cyclic(origin_segment)) {
      if (dsl.size() % 2 == 0) sign_ratio *= -1;
    } else if (is_cyclic(dsl.back()))
      sign_ratio *= -1;
    dsl.insert(destination_it, origin_segment);

    // Remove the segment at origin and update sign
    auto &sl = config.seglists[origin_color];
    if (is_cyclic(origin_segment)) {
      if (sl.size() % 2 == 1) sign_ratio *= -1;
    } else if (is_cyclic(sl.back()))
      sign_ratio *= -1;
    sl.erase(sl.begin() + origin_index);

    return sign_ratio;
  }

  //--------------------------------------------------
  void move_segment::reject() {
    SPDLOG_TRACE("\n - - - - - ====> REJECT - - - - - - - - - - -\n");
    wdata.dets[origin_color].reject_last_try();
    wdata.dets[destination_color].reject_last_try();
  }
}; // namespace moves
