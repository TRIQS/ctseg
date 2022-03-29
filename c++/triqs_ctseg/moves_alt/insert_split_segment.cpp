#include "insert_split_segment.hpp"
#include "../logs.hpp"
#include <cmath>

namespace moves {

  // Find the segement to the left of tau (accounting for cyclicity), and whether tau is inside
  void insert_split_segment::prep_insertion(std::vector<segment_t> const &seglist, dimtime_t tau) {
    is_inside = false;
    seg_it    = seglist.cend();
    if (seglist.empty()) return;
    seg_it = find_segment_left(seglist, segment_t{tau, wdata.qmc_zero});
    if (seg_it->tau_c > tau) { // Actually found segment to the left
      if (seg_it->tau_cdag == tau or seg_it->tau_c == tau) throw;
      if (seg_it->tau_cdag < tau or is_cyclic(*seg_it)) is_inside = true;
    }
    // Found segment is in fact to the right
    else {
      if (seg_it->tau_c == tau or seglist.back().tau_cdag == tau) throw;
      seg_it = --seglist.cend();
      if (is_cyclic(seglist.back()) and seg_it->tau_cdag < tau) is_inside = true;
    }
  }

  // ------------------------------------------------------------------

  double insert_split_segment::attempt() {

#if 0
    LOG("\n =================== ATTEMPT INSERT/SPLIT ================ \n");

    // Select insertion color
    color    = rng(wdata.n_color);
    auto &sl = config.seglists[color];
    LOG("Inserting/splitting at color {}", color);

    // Select first time
    tau1 = dimtime_t::random(rng, wdata.qmc_zero, wdata.qmc_beta);
    LOG("First time is {}", tau1);
    try {
      prep_insertion(sl, tau1);
    } catch (...) {
      LOG("First time coincides with an operator.");
      return 0;
    }

    seg_idx            = std::distance(sl.cbegin(), seg_it);
    double trace_ratio = 0, det_ratio = 0, prop_ratio = 0, window_length = 0;

    // ================== SPLIT ===================
    if (is_inside) {
      LOG("First time is inside a segment: attempt SPLIT at position {}.", seg_idx);

      // ------------- Choose second time ---------
      tau2 = dimtime_t::random(rng, seg_it->tau_cdag, seg_it->tau_c);
      if (tau2 == tau1) {
        LOG("Second time equal to first time.");
        return 0;
      }
      splitting_full_line = is_full_line(*seg_it, wdata.fac);
      if (not splitting_full_line and tau1 - seg_it->tau_c > tau2 - seg_it->tau_c)
        std::swap(tau1, tau2); // Unless splitting full line, tau1 is left of tau2.
      LOG("Split: adding c at {}, cdag at {}", tau2, tau1);
      window_length        = seg_it->length();
      auto removed_segment = segment_t{tau1, tau2}; // "antisegment" : careful with order of c, cdag

      // ------------  Trace ratio  -------------
      double ln_trace_ratio = 0;
      for (auto c : range(wdata.n_color)) {
        if (c != color) {
          ln_trace_ratio -= -wdata.U(color, c) * overlap(config.seglists[c], removed_segment, wdata.fac);
          ln_trace_ratio -= wdata.mu(c) * removed_segment.length();
          if (wdata.has_Dt)
            ln_trace_ratio -= K_overlap(config.seglists[c], removed_segment, slice_target_to_scalar(wdata.K, color, c));
        }
      }
      trace_ratio = std::exp(ln_trace_ratio);

      // ------------  Det ratio  ---------------

      /* We insert tau_cdag as a line (first index) and tau_c as a column (second index). The index always corresponds to the
     segment the tau_c/tau_cdag belongs to. Here, the cdag is always inserted at the position of the segment we are splitting.
     The c insertion position depends on whether we are splitting a full line and whether the cut out segment is cyclic, computed 
     in right_seg_idx. */

      right_seg_idx = (splitting_full_line or is_cyclic(removed_segment)) ? 0 : seg_idx + 1;
      det_ratio     = wdata.dets[color].try_insert(seg_idx, right_seg_idx, {tau1, 0}, {tau2, 0});

    } // Split

    // ==================== INSERT =====================
    else {
      LOG("First time is not in a segment: attempt INSERT to the right of position {}.", seg_idx);

      // ------------ Find insertion window [wtau1,wtau2] ------------
      insert_into_empty_line = seg_it == sl.end();
      dimtime_t wtau1, wtau2;
      if (insert_into_empty_line) {
        wtau1         = wdata.qmc_beta;
        wtau2         = wdata.qmc_zero;
        window_length = wdata.beta;
      } else {
        wtau1            = seg_it->tau_cdag;
        auto next_seg_it = (seg_it == --sl.end()) ? sl.begin() : seg_it++;
        wtau2            = next_seg_it->tau_c;
        window_length    = double(wtau1 - wtau2);
      }

      // ------------- Choose second time --------------
      tau2 = dimtime_t::random(rng, wtau2, wtau1);
      if (tau2 == tau1) {
        LOG("Second time equal to first time.");
        return 0;
      }
      // Unless inserting into empty line, tau1 is left of tau2.
      if (not insert_into_empty_line and tau1 - seg_it->tau_cdag > tau2 - seg_it->tau_cdag) std::swap(tau1, tau2);

      proposed_segment             = segment_t{tau1, tau2};
      auto proposed_segment_length = double(tau1 - tau2); // can be cyclic.
      insert_it                    = ++seg_it;
      if (seg_it == --sl.end() and tau1 > seg_it->tau_cdag) insert_it = sl.cbegin();

      LOG("Inserting c at {}, cdag at {}", tau1, tau2);

      // ------------  Trace ratio  -------------
      double ln_trace_ratio = 0;
      for (auto c : range(wdata.n_color)) {
        if (c != color) {
          ln_trace_ratio += -wdata.U(color, c) * overlap(config.seglists[c], proposed_segment, wdata.fac);
          ln_trace_ratio += wdata.mu(c) * proposed_segment_length;
          if (wdata.has_Dt)
            ln_trace_ratio +=
               K_overlap(config.seglists[c], proposed_segment, slice_target_to_scalar(wdata.K, color, c));
        }
      }
      trace_ratio = std::exp(ln_trace_ratio);

      // ------------  Det ratio  -------------------
      // insert_idx is the position of the proposed segment if inserted
      insert_idx = std::distance(sl.cbegin(), insert_it);
      // We insert tau_cdag as a line (first index) and tau_c as a column (second index). The index always corresponds to the
      // segment the tau_c/tau_cdag belongs to.
      det_ratio = wdata.dets[color].try_insert(insert_idx, insert_idx, {proposed_segment.tau_cdag, 0},
                                               {proposed_segment.tau_c, 0});

    } // Insert

    // -------------- Prop ratio --------------------
    bool no_swap                  = splitting_full_line or insert_into_empty_line;
    double future_number_segments = splitting_full_line ? 1 : sl.size() + 1.0;
    prop_ratio                    = wdata.beta * window_length / ((no_swap ? 1 : 2) * future_number_segments);

    LOG("trace_ratio  = {}, prop_ratio = {}, det_ratio = {}", trace_ratio, prop_ratio, det_ratio);

    double prod = trace_ratio * det_ratio * prop_ratio;
    det_sign    = (det_ratio > 0) ? 1.0 : -1.0;
    return (std::isfinite(prod) ? prod : det_sign);
#endif
    return 0;
  }

  //--------------------------------------------------

  double insert_split_segment::accept() {

    LOG("\n - - - - - ====> ACCEPT - - - - - - - - - - -\n");

    // Insert the times into the det
    wdata.dets[color].complete_operation();

    double sign_ratio = 1;
    auto &sl          = config.seglists[color];

    if (is_inside) { // Split
      // Split the segment and compute the sign ratio
      if (splitting_full_line) {
        auto new_segment = segment_t{tau2, tau1};
        if (is_cyclic(new_segment)) sign_ratio = -1;
        sl[seg_idx] = new_segment;
      } else {
        auto new_segment_left  = segment_t{seg_it->tau_c, tau1};
        auto new_segment_right = segment_t{tau2, seg_it->tau_cdag};
        if (is_cyclic(*seg_it)) {
          bool kept_number_cyclic = is_cyclic(new_segment_left) or is_cyclic(new_segment_right);
          // If we are splitting a cyclic segment and both new segments are not cyclic, get a - sign
          if (not kept_number_cyclic) sign_ratio = -1;
        }
        // Update the proposed segment
        sl[seg_idx] = new_segment_left;
        // Insert a new segment
        sl.insert(sl.begin() + right_seg_idx, new_segment_right);
      }
    } else { // Insert
      // Compute the sign ratio
      if (is_cyclic(proposed_segment)) sign_ratio = -1;
      // Insert the segment in an ordered list
      sl.insert(insert_it, proposed_segment);
    }

    LOG("Sign ratio is {}", sign_ratio);

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
  void insert_split_segment::reject() {
    LOG("\n - - - - - ====> REJECT - - - - - - - - - - -\n");
    wdata.dets[color].reject_last_try();
  }
}; // namespace moves
