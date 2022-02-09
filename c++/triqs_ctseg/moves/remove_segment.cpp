#include "remove_segment.hpp"
namespace moves {

  // FIXME : LIGNE PLEINE ?? dans overlap ???

  double remove_segment::attempt() {

    SPDLOG_LOGGER_TRACE("\n =================== ATTEMPT REMOVE ================ \n");

    // ------------ Choice of segment --------------

    // Select removal color
    color = rng(data.n_color);
    auto &sl       = config.seglist[color];
    SPDLOG_LOGGER_TRACE("Removing at color {}", n_color);

    // If color is empty, nothing to remove
    if sl.empty() return 0; 

    // Select segment to remove 
    proposed_segment_index = rng(sl.size());
    proposed_segment = *proposed_segment_index;

    SPDLOG_LOGGER_TRACE("Removing c at{}, cdag at {}", proposed_segment.tau_c, proposed_segment.tau_cdag);

    // ------------  Trace ratio  -------------
        // FIXME : here we will need the K function integral 
    double ln_trace_ratio = 0;
    for (int c : seglists) {
      if (c != color) ln_trace_ratio -= overlap(proposed_segment, config.seglists[c]);
    }
    double trace_ratio = std::exp(ln_trace_ratio);

    // ------------  Det ratio  ---------------

    // FIXME

    // ------------  Proposition ratio ------------

    double current_number_segments = sl.size();
    double future_number_segments = std::max(sl.size()-1,1); 
    // Limits of insertion interval for reverse move, initialise at (beta,0)
    qmc_time_t tau_left = time_point_factory.get_upper_pt();
    qmc_time_t tau_right = time_point_factory.get_lower_pt();
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

    SPDLOG_LOGGER_TRACE("\n - - - - - ====> ACCEPT - - - - - - - - - - -\n");

    data.dets[color].complete_operation();
    // Remove the segment
    config.seglists[color].erase(proposed_segment_insert_pos);

    // FIXME ??? SIGNE ???
    double sign_ratio = 1; 
    return sign_ratio;
  }

  //--------------------------------------------------
  void remove_segment::reject() {
    SPDLOG_LOGGER_TRACE("\n - - - - - ====> REJECT - - - - - - - - - - -\n");
    data.dets[color].reject_last_try();
  }
};
} // namespace moves
