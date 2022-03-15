#include "regroup_segment.hpp"
#include "../logs.hpp"

namespace moves {

  double regroup_segment::attempt() {

    SPDLOG_TRACE("\n =================== ATTEMPT REGROUP ================ \n");

    // ------------ Choice of segment --------------
    // Select color
    color    = rng(wdata.n_color);
    auto &sl = config.seglists[color];
    SPDLOG_TRACE("Regrouping at color {}", color);

    // If no segments nothing to regroup
    if (sl.empty()) return 0;

    // Select pair of segments (or cyclic segment) to regroup
    making_full_line = sl.size() == 1;
    if (making_full_line) {
      if (is_full_line(sl[0], fac)) return 0; // If segment is full line nothing to regroup
      left_segment_index  = 0;
      right_segment_index = 0;
    } else {
      left_segment_index  = rng(sl.size());
      right_segment_index = (left_segment_index == sl.size() - 1) ? 0 : left_segment_index + 1;
    }
    left_segment  = sl[left_segment_index];
    right_segment = sl[right_segment_index];

    auto inserted_segment =
       segment_t{left_segment.tau_cdag, right_segment.tau_c}; // "antisegment" : careful with order of c, cdag
    auto inserted_segment_length = double(inserted_segment.tau_c - inserted_segment.tau_cdag);

    SPDLOG_TRACE("Regroup: removing c at {}, cdag at {}", right_segment.tau_c, left_segment.tau_cdag);

    // ------------  Trace ratio  -------------
    double ln_trace_ratio = 0;
    for (auto c : range(wdata.n_color)) {
      if (c != color) {
        ln_trace_ratio += -wdata.U(color, c) * overlap(config.seglists[c], inserted_segment, fac);
        ln_trace_ratio += wdata.mu(c) * inserted_segment_length;
        if (wdata.has_Dt)
          ln_trace_ratio += K_overlap(config.seglists[c], inserted_segment, slice_target_to_scalar(wdata.K, color, c));
      }
    }
    double trace_ratio = std::exp(ln_trace_ratio);

    // ------------  Det ratio  ---------------
    // We remove a cdag (first index) from the left segment and a c (second index) from the right segment.
    auto det_ratio = wdata.dets[color].try_remove(
       left_segment_index, right_segment_index); // FIXME: ordering in det when regrouping into cyclic segment

    // additional sign due to "roll" is not computed at this stage, cf accept

    // ------------  Proposition ratio ------------

    double future_number_segments   = making_full_line ? 1 : int(sl.size()) - 1;
    double current_number_intervals = sl.size();
    // Length of future segment
    qmc_time_t l = fac.get_upper_pt();
    if (not making_full_line) l = left_segment.tau_c - right_segment.tau_cdag;
    double prop_ratio = (future_number_segments * l * l / (making_full_line ? 1 : 2)) / current_number_intervals;

    SPDLOG_TRACE("trace_ratio  = {}, prop_ratio = {}, det_ratio = {}", trace_ratio, prop_ratio, det_ratio);

    return trace_ratio * det_ratio * prop_ratio;
  }

  //--------------------------------------------------

  double regroup_segment::accept() {

    SPDLOG_TRACE("\n - - - - - ====> ACCEPT - - - - - - - - - - -\n");

    // Update the dets
    wdata.dets[color].complete_operation();

    // Regroup segments and copute the sign ratio
    auto &sl          = config.seglists[color];
    double sign_ratio = 1;
    if (making_full_line) {
      if (is_cyclic(sl[left_segment_index])) sign_ratio = -1;
      sl[left_segment_index] = segment_t{wdata.qmc_beta, wdata.qmc_zero};
    } else {
      auto new_segment        = segment_t{left_segment.tau_c, right_segment.tau_cdag};
      bool old_segment_cyclic = is_cyclic(sl[left_segment_index]) or is_cyclic(sl[right_segment_index]);
      // If one of the regrouped segments was cyclic, we still have a cyclic segment and we changed
      // the parity of the number of segments
      if (old_segment_cyclic) sign_ratio = -1;
      // Otherwise, we may have obtained a cyclic segment from two non-cyclic segments
      else if (is_cyclic(new_segment)) {
        // In this case, a cdag that belonged to the first segment now belongs to the last segment:
        // we need to "roll up" the det and this produces an additional sign.
        sign_ratio *= wdata.dets[color].roll_matrix(work_data_t::det_t::Up);
        // A further additional sign if we got to an odd configuration.
        if (sl.size() % 2 == 0) sign_ratio *= -1;
      }
      // Update the left segment
      sl[left_segment_index] = new_segment;
      // Remove the "other" segment
      sl.erase(sl.begin() + right_segment_index);
    }
    return sign_ratio;
  }

  //--------------------------------------------------
  void regroup_segment::reject() {
    SPDLOG_TRACE("\n - - - - - ====> REJECT - - - - - - - - - - -\n");
    wdata.dets[color].reject_last_try();
  }
} // namespace moves
