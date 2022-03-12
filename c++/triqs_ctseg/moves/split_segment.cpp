#include "split_segment.hpp"

namespace moves {

  double split_segment::attempt() {

    SPDLOG_LOGGER_TRACE("\n =================== ATTEMPT SPLIT ================ \n", void);

    // ------------ Choice of segment --------------
    // Select color
    color    = rng(wdata.n_color);
    auto &sl = config.seglists[color];
    SPDLOG_LOGGER_TRACE("Splitting at color {}", color);

    // If color is empty, nothing to split
    if (sl.empty()) return 0;

    // Select segment to split
    proposed_segment_idx = rng(sl.size());
    proposed_segment     = sl[proposed_segment_idx];
    // Select splitting points (tau_left,tau_right)
    auto qmc_zero  = fac.get_lower_pt();
    qmc_time_t l   = proposed_segment.tau_c - proposed_segment.tau_cdag;
    qmc_time_t dt1 = fac.get_random_pt(rng, qmc_zero, l);
    qmc_time_t dt2 = fac.get_random_pt(rng, qmc_zero, l);
    if (dt1 == dt2) return 0;
    if (dt1 == qmc_zero || dt2 == qmc_zero || dt1 == l || dt2 == l) return 0;
    full_line = is_full_line(proposed_segment, fac);
    if (dt1 > dt2 && !full_line)
      std::swap(dt1, dt2); // If splitting a full line, the order of tau_left and tau_right is not fixed
    tau_left                    = proposed_segment.tau_c - dt1; // dt1 < dt2
    tau_right                   = proposed_segment.tau_c - dt2;
    auto removed_segment        = segment_t{tau_left, tau_right}; // "antisegment" : careful with order of c, cdag
    auto removed_segment_length = double(tau_left - tau_right);   // accounts for cyclicity
    bool removed_segment_cyclic = tau_left < tau_right;
    right_segment_idx           = (full_line or removed_segment_cyclic) ? 0 : proposed_segment_idx + 1;

    SPDLOG_LOGGER_TRACE("Split: adding c at {}, cdag at {}", tau_right, tau_left);

    // ------------  Trace ratio  -------------
    double ln_trace_ratio = 0;
    for (auto c : range(wdata.n_color)) {
      if (c != color) {
        ln_trace_ratio -= -wdata.U(color, c) * overlap(config.seglists[c], removed_segment, fac);
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

    SPDLOG_LOGGER_TRACE("trace_ratio  = {}, prop_ratio = {}, det_ratio = {}", trace_ratio, prop_ratio, det_ratio);

    return trace_ratio * det_ratio * prop_ratio;
  }

  //--------------------------------------------------

  double split_segment::accept() {

    SPDLOG_LOGGER_TRACE("\n - - - - - ====> ACCEPT - - - - - - - - - - -\n", void);

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
        // If we are splitting a cyclic segment and one of the new segments is cyclic, we changed the parity
        // of the number of segments and get a - sign
        if (kept_number_cyclic) sign_ratio = -1;
        // Otherwise, we get a - sign only if we were in an odd configuration
        else if (sl.size() % 2 == 1)
          sign_ratio = -1;
        // Always get a - sign if the config had a cyclic segment and we did not touch it
      } else if (is_cyclic(sl.back()))
        sign_ratio = -1;
      // Update the proposed segment
      sl[proposed_segment_idx] = new_segment_left;
      // Insert a new segment
      sl.insert(sl.begin() + right_segment_idx, new_segment_right);
    }

    return sign_ratio;
  }

  //--------------------------------------------------
  void split_segment::reject() {
    SPDLOG_LOGGER_TRACE("\n - - - - - ====> REJECT - - - - - - - - - - -\n", void);
    wdata.dets[color].reject_last_try();
  }
}; // namespace moves
