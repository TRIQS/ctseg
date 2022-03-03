#include "remove_segment.hpp"
namespace moves {

  double remove_segment::attempt() {

    SPDLOG_LOGGER_TRACE("\n =================== ATTEMPT REMOVE ================ \n",void);

    // ------------ Choice of segment --------------

    // Select removal color
    color = rng(wdata.n_color);
    auto &sl       = config.seglists[color];
    SPDLOG_LOGGER_TRACE("Removing at color {}", color);

    // If color is empty, nothing to remove
    if (sl.empty()) return 0; 

    // Select segment to remove 
    proposed_segment_index = rng(sl.size());
    proposed_segment = sl[proposed_segment_index];
    auto qmc_beta = time_point_factory.get_upper_pt();
    auto qmc_zero = time_point_factory.get_lower_pt();
    if (proposed_segment.tau_c == qmc_beta && proposed_segment.tau_cdag == qmc_zero) return 0; // If segment is a full line do not remove 
    auto proposed_segment_length = double(proposed_segment.tau_c-proposed_segment.tau_cdag);

    SPDLOG_LOGGER_TRACE("Removing c at{}, cdag at {}", proposed_segment.tau_c, proposed_segment.tau_cdag);

    // ------------  Trace ratio  -------------
    double ln_trace_ratio = 0;
    for (auto c : range(wdata.n_color)) {
      if (c != color) {
        ln_trace_ratio -= -wdata.U(color,c)*overlap(config.seglists[c],proposed_segment,time_point_factory);
        ln_trace_ratio -= wdata.mu(c)*proposed_segment_length;
        if (wdata.has_Dt) ln_trace_ratio -= K_overlap(config.seglists[c],proposed_segment,slice_target_to_scalar(wdata.K,color,c)); // FIXME. Is the syntax right for slice????
      }
    }
    double trace_ratio = std::exp(ln_trace_ratio);

    // ------------  Det ratio  ---------------

    // FIXME
    double det_ratio = 0;

    // ------------  Proposition ratio ------------

    double current_number_segments = sl.size();
    double future_number_segments = std::max(int(sl.size()-1),1); 
    // Limits of insertion interval for reverse move, initialise at (beta,0)
    qmc_time_t tau_left = qmc_beta;
    qmc_time_t tau_right = qmc_zero;
    if (current_number_segments != 1) {
        bool is_last_segment = proposed_segment_index == sl.size() - 1; 
        bool is_first_segment = proposed_segment_index == 0; 
        // Look at segments on either side of proposed_segment, accounting for cyclicity 
        tau_right = sl[is_last_segment ? 0 : proposed_segment_index + 1].tau_c; 
        tau_left = sl[is_first_segment ? sl.size() - 1 : proposed_segment_index - 1].tau_cdag; 
    }
    qmc_time_t l = tau_left - tau_right; 
    double prop_ratio = (future_number_segments * l * l) / current_number_segments;

    SPDLOG_LOGGER_TRACE("trace_ratio  = {}, prop_ratio = {}, det_ratio = {}", trace_ratio, prop_ratio, det_ratio);

    return trace_ratio * det_ratio * prop_ratio;
  }

  //--------------------------------------------------

  double remove_segment::accept() {

    SPDLOG_LOGGER_TRACE("\n - - - - - ====> ACCEPT - - - - - - - - - - -\n",void);

    //data.dets[color].complete_operation();
    // Remove the segment
    auto &sl = config.seglists[color];
    auto proposed_segment_it = std::next(sl.begin(),proposed_segment_index);
    sl.erase(proposed_segment_it);

    // FIXME ??? SIGNE ???
    double sign_ratio = 1; 
    return sign_ratio;
  }

  //--------------------------------------------------
  void remove_segment::reject() {
    SPDLOG_LOGGER_TRACE("\n - - - - - ====> REJECT - - - - - - - - - - -\n",void);
    //data.dets[color].reject_last_try();
  }
} // namespace moves
