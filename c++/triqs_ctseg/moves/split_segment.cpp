#include "split_segment.hpp"
#include "../logs.hpp"

namespace moves {

  double split_segment::attempt() {

    LOG("\n =================== ATTEMPT SPLIT ================ \n");

    // ------------ Choice of segment --------------
    // Select color
    color    = rng(wdata.n_color);
    auto &sl = config.seglists[color];
    LOG("Splitting at color {}", color);

    // If color is empty, nothing to split
    if (sl.empty()) {
      LOG("Line is empty");
      return 0;
    }
    // Select segment to split
    proposed_segment_idx = rng(sl.size());
    proposed_segment     = sl[proposed_segment_idx];
    full_line            = is_full_line(proposed_segment, fac);
    if (full_line) LOG("Splitting full line.");
    // Select splitting points (tau_left,tau_right)
    qmc_time_t l   = full_line ? wdata.qmc_beta : proposed_segment.tau_c - proposed_segment.tau_cdag;
    qmc_time_t dt1 = fac.get_random_pt(rng, wdata.qmc_zero, l);
    qmc_time_t dt2 = fac.get_random_pt(rng, wdata.qmc_zero, l);
    if (dt1 == dt2 || dt1 == wdata.qmc_zero || dt2 == wdata.qmc_zero || dt1 == l || dt2 == l) {
      LOG("Generated time at boundary");
      return 0;
    }
    if (dt1 > dt2 && !full_line)
      std::swap(dt1, dt2); // If splitting a full line, the order of tau_left and tau_right is not fixed
    tau_left                    = proposed_segment.tau_c - dt1; // dt1 < dt2
    tau_right                   = proposed_segment.tau_c - dt2;
    auto removed_segment        = segment_t{tau_left, tau_right}; // "antisegment" : careful with order of c, cdag
    auto removed_segment_length = double(tau_left - tau_right);   // accounts for cyclicity
    bool removed_segment_cyclic = tau_left < tau_right;
    right_segment_idx           = (full_line or removed_segment_cyclic) ? 0 : proposed_segment_idx + 1;

    LOG("Split: adding c at {}, cdag at {}", tau_right, tau_left);

    // ------------  Trace ratio  -------------
    double ln_trace_ratio = 0;
    for (auto c : range(wdata.n_color)) {
      if (c != color) {
        ln_trace_ratio -= -wdata.U(color, c) * overlap(config.seglists[c], removed_segment, fac);
        LOG("Overlap is {}", ln_trace_ratio);
        ln_trace_ratio -= wdata.mu(c) * removed_segment_length;
        if (wdata.has_Dt)
          ln_trace_ratio -= K_overlap(config.seglists[c], removed_segment, slice_target_to_scalar(wdata.K, color, c));
      }
    }
    double trace_ratio = std::exp(ln_trace_ratio);

    // ------------  Det ratio  ---------------

    /* We insert tau_cdag as a line (first index) and tau_c as a column (second index). The index always corresponds to the
     segment the tau_c/tau_cdag belongs to. Here, the cdag is always inserted at the position of the segment we are splitting.
     The c insertion position depends on whether we are splitting a full line and whether the cut out segment is cyclic, computed 
     in right_segment_idx. */
    auto det_ratio =
       wdata.dets[color].try_insert(proposed_segment_idx, right_segment_idx, {tau_left, 0}, {tau_right, 0});

    // ------------  Proposition ratio ------------

    double current_number_segments = sl.size();
    double future_number_intervals = full_line ? 1 : double(sl.size()) + 1.0;
    double prop_ratio = (future_number_intervals) / (current_number_segments * l * l / (full_line ? 1 : 2));

    LOG("trace_ratio  = {}, prop_ratio = {}, det_ratio = {}", trace_ratio, prop_ratio, det_ratio);

    double prod = trace_ratio * det_ratio * prop_ratio;
    det_sign    = (det_ratio > 0) ? 1.0 : -1.0;
    return (std::isfinite(prod) ? prod : det_sign);
  }

  //--------------------------------------------------

  double split_segment::accept() {

    LOG("\n - - - - - ====> ACCEPT - - - - - - - - - - -\n");

    // Update the dets
    wdata.dets[color].complete_operation();

    auto &sl          = config.seglists[color];
    double sign_ratio = 1;
    // Split the segment and compute the sign ratio
    if (full_line) {
      auto new_segment = segment_t{tau_right, tau_left};
      if (is_cyclic(new_segment)) sign_ratio = -1;
      sl[proposed_segment_idx] = new_segment;
    } else {
      auto new_segment_left  = segment_t{proposed_segment.tau_c, tau_left};
      auto new_segment_right = segment_t{tau_right, proposed_segment.tau_cdag};
      if (is_cyclic(proposed_segment)) {
        bool kept_number_cyclic = is_cyclic(new_segment_left) or is_cyclic(new_segment_right);
        // If we are splitting a cyclic segment and both new segments are not cyclic, get a - sign
        if (not kept_number_cyclic) sign_ratio = -1;
      }
      LOG("Sign ratio is {}", sign_ratio);

      // Update the proposed segment
      sl[proposed_segment_idx] = new_segment_left;
      // Insert a new segment
      sl.insert(sl.begin() + right_segment_idx, new_segment_right);
    }

    // Check invariant
    ALWAYS_EXPECTS((sign_ratio * det_sign == 1.0),
                   "Error: move has produced negative sign! Det sign is {} and additional sign is {}.", det_sign,
                   sign_ratio);
#ifdef CHECK_INVARIANTS
    check_invariant(config, wdata.dets);
#endif
    SPDLOG_TRACE("Configuration is {}", config);

    return sign_ratio;
  }

  //--------------------------------------------------
  void split_segment::reject() {
    LOG("\n - - - - - ====> REJECT - - - - - - - - - - -\n");
    wdata.dets[color].reject_last_try();
  }
}; // namespace moves
