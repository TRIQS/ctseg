#include "insert_segment.hpp"
#include "../logs.hpp"
#include <cmath>

namespace moves {

  double insert_segment::attempt() {

    LOG("\n =================== ATTEMPT INSERT ================ \n");

    // ------------ Choice of segment --------------

    // Select insertion color
    color    = rng(wdata.n_color);
    auto &sl = config.seglists[color];
    LOG("Inserting at color {}", color);

    // Select insertion window [wtau_left,wtau_right]
    qmc_time_t wtau_left = wdata.qmc_beta, wtau_right = wdata.qmc_zero;

    if (not sl.empty()) {
      if (is_full_line(sl.back(), fac)) {
        LOG("Full line, cannot insert.");
        return 0;
      }
      // Randomly choose one existing segment
      long seg_idx         = rng(sl.size());
      wtau_left            = sl[seg_idx].tau_cdag; // wtau_left is cdag of this segment
      bool is_last_segment = seg_idx == sl.size() - 1;
      wtau_right = sl[is_last_segment ? 0 : seg_idx + 1].tau_c; // wtau_right is c of next segment, possibly cyclic
    }

    LOG("Insertion window is wtau_left = {}, wtau_right = {}", wtau_left, wtau_right);
    qmc_time_t window_length = sl.empty() ? wdata.qmc_beta : wtau_left - wtau_right;

    // Choose two random times in insertion window
    auto dt1 = wdata.make_random_time(rng, window_length);
    auto dt2 = wdata.make_random_time(rng, window_length);
    if (dt1 == dt2) {
      LOG("Generated equal times");
      return 0;
    }
    if (dt1 > dt2 and !sl.empty()) std::swap(dt1, dt2); // if inserting into an empty line, two ways to insert
    prop_seg = segment_t{wtau_left - dt1, wtau_left - dt2};
    // The index of prop_seg if it is inserted in the list of segments.
    prop_seg_it       = std::upper_bound(sl.begin(), sl.end(), prop_seg);
    long prop_seg_idx = std::distance(sl.begin(), prop_seg_it);

    LOG("Inserting segment at position {}, with c at {}, cdag at {}", prop_seg_idx, prop_seg.tau_c, prop_seg.tau_cdag);

    // ------------  Trace ratio  -------------
    double ln_trace_ratio = wdata.mu(color) * prop_seg.length();
    for (auto c : range(wdata.n_color)) {
      if (c != color) {
        ln_trace_ratio += -wdata.U(color, c) * overlap(config.seglists[c], prop_seg, fac);
        if (wdata.has_Dt)
          ln_trace_ratio += K_overlap(config.seglists[c], prop_seg, slice_target_to_scalar(wdata.K, color, c));
      }
    }
    double trace_ratio = std::exp(ln_trace_ratio);

    // ------------  Det ratio  ---------------
    // We insert tau_cdag as a line (first index) and tau_c as a column (second index). The index always corresponds to the
    // segment the tau_c/tau_cdag belongs to.
    auto det_ratio =
       wdata.dets[color].try_insert(prop_seg_idx, prop_seg_idx, {prop_seg.tau_cdag, 0}, {prop_seg.tau_c, 0});

    // ------------  Proposition ratio ------------

    double current_number_intervals = std::max(1.0, double(sl.size()));
    double future_number_segments   = double(sl.size()) + 1;
    double prop_ratio =
       (current_number_intervals * window_length * window_length / (sl.empty() ? 1 : 2)) / future_number_segments;
    // Account for absence of time swapping when inserting into empty line.

    LOG("trace_ratio  = {}, prop_ratio = {}, det_ratio = {}", trace_ratio, prop_ratio, det_ratio);

    double prod = trace_ratio * det_ratio * prop_ratio;
    det_sign    = (det_ratio > 0) ? 1.0 : -1.0;
    return (std::isfinite(prod) ? prod : det_sign);
  }

  //--------------------------------------------------

  double insert_segment::accept() {

    LOG("\n - - - - - ====> ACCEPT - - - - - - - - - - -\n");

    // Insert the times into the det
    wdata.dets[color].complete_operation();

    // Compute the sign ratio
    double sign_ratio = 1;
    auto &sl          = config.seglists[color];
    if (is_cyclic(prop_seg)) sign_ratio = -1.0;
    LOG("Sign ratio is {}", sign_ratio);

    // Insert the segment in an ordered list
    sl.insert(prop_seg_it, prop_seg);

    // Check invariant
#ifdef CHECK_INVARIANTS
    check_invariant(config, wdata.dets);
#endif
    ALWAYS_EXPECTS((sign_ratio * det_sign == 1.0),
                   "Error: move has produced negative sign! Det sign is {} and additional sign is {}.", det_sign,
                   sign_ratio);
    SPDLOG_TRACE("Configuration is {}", config);

    return sign_ratio;
  }

  //--------------------------------------------------
  void insert_segment::reject() {
    LOG("\n - - - - - ====> REJECT - - - - - - - - - - -\n");
    wdata.dets[color].reject_last_try();
  }
}; // namespace moves
