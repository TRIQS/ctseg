#include "split_segment.hpp"

namespace moves {

  double split_segment::attempt() {

    SPDLOG_LOGGER_TRACE("\n =================== ATTEMPT SPLIT ================ \n",void);

    // ------------ Choice of segment --------------
    // Select color
    color = rng(wdata.n_color);
    auto &sl       = config.seglists[color];
    SPDLOG_LOGGER_TRACE("Splitting at color {}", color);

    // If color is empty, nothing to split
    if (sl.empty()) return 0; 

    // Select segment to split 
    proposed_segment_index = rng(sl.size());
    proposed_segment = sl[proposed_segment_index];
    // Select splitting points (tau_left,tau_right)
    auto qmc_zero = time_point_factory.get_lower_pt();
    qmc_time_t l = proposed_segment.tau_c - proposed_segment.tau_cdag; 
    qmc_time_t dt1 = time_point_factory.get_random_pt(rng,qmc_zero,l); 
    qmc_time_t dt2 = time_point_factory.get_random_pt(rng,qmc_zero,l);
    if (dt1 == dt2) return 0; 
    if (dt1 == qmc_zero || dt2 == qmc_zero || dt1 == l || dt2 == l) return 0; 
    full_line = is_full_line(proposed_segment,time_point_factory);
    if (dt1 > dt2 && !full_line) std::swap(dt1,dt2); // If splitting a full line, the order of tau_left and tau_right is not fixed
    tau_left = proposed_segment.tau_c - dt1; 
    tau_right = proposed_segment.tau_c - dt2; 
    auto removed_segment = segment_t{tau_left,tau_right};
    auto removed_segment_length = double(tau_left - tau_right);

    SPDLOG_LOGGER_TRACE("Split: adding c at {}, cdag at {}", tau_right, tau_left);

    // ------------  Trace ratio  -------------
    double ln_trace_ratio = 0;
    for (auto c : range(wdata.n_color)) {
      if (c != color) {
        ln_trace_ratio -= -wdata.U(color,c)*overlap(config.seglists[c],removed_segment,time_point_factory);
        ln_trace_ratio -= wdata.mu(c)*removed_segment_length; 
        if (wdata.has_Dt) ln_trace_ratio -= K_overlap(config.seglists[c],removed_segment,slice_target_to_scalar(wdata.K,color,c)); // FIXME. Is the syntax right for slice????
      }
    }
    double trace_ratio = std::exp(ln_trace_ratio);

    // ------------  Det ratio  ---------------

    // FIXME
    double det_ratio = 0;

    // ------------  Proposition ratio ------------

    double current_number_segments = full_line ? 2 : int(sl.size()); // Account for the two ways of splitting a full line 
    double future_number_intervals = full_line ? 1 : int(sl.size()) + 1; 
    double prop_ratio = (future_number_intervals) / (current_number_segments * l * l);

    SPDLOG_LOGGER_TRACE("trace_ratio  = {}, prop_ratio = {}, det_ratio = {}", trace_ratio, prop_ratio, det_ratio);

    return trace_ratio * det_ratio * prop_ratio;
  }

  //--------------------------------------------------

  double split_segment::accept() {

    SPDLOG_LOGGER_TRACE("\n - - - - - ====> ACCEPT - - - - - - - - - - -\n",void);

    //data.dets[color].complete_operation();
    // Split the segment
    auto &sl = config.seglists[color];
    if (is_full_line(proposed_segment,time_point_factory)) {
      auto new_segment = segment_t{tau_right,tau_left};
      sl[proposed_segment_index] = new_segment; 
    }
    else {
      auto new_segment_left = segment_t{proposed_segment.tau_c,tau_left}; 
      auto new_segment_right = segment_t{tau_right,proposed_segment.tau_cdag}; 
      auto proposed_segment_it = std::next(sl.begin(),proposed_segment_index);
      sl.insert(proposed_segment_it,new_segment_left); 
      sl[proposed_segment_index + 1] = new_segment_right; 
    }

    // FIXME ??? SIGNE ???
    double sign_ratio = 1; 
    return sign_ratio;
  }

  //--------------------------------------------------
  void split_segment::reject() {
    SPDLOG_LOGGER_TRACE("\n - - - - - ====> REJECT - - - - - - - - - - -\n",void);
    //data.dets[color].reject_last_try();
  }
}; // namespace moves
