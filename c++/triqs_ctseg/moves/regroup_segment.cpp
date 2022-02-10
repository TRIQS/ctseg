#include "regroup_segment.hpp"

namespace moves {

  double regroup_segment::attempt() {

    SPDLOG_LOGGER_TRACE("\n =================== ATTEMPT SPLIT ================ \n");

    // ------------ Choice of segment --------------
    // Select color
    color = rng(data.n_color);
    auto &sl       = config.seglist[color];
    SPDLOG_LOGGER_TRACE("Regrouping at color {}", n_color);

    // If color has 0 or 1 segment; nothing to regroup
    if (sl.empty() || sl.size() == 1) return 0; 

    // Select pair of segments to regroup 
    left_segment_index = rng(sl.size());
    right_segment_index = (left_segment_index == sl.size() - 1) ? 0 : left_segment_index + 1; 
    left_segment = sl[left_segment_index];
    right_segment = sl[right_segment_index];

    SPDLOG_LOGGER_TRACE("Regroup: removing c at {}, cdag at {}", right_segment.tau_c, left_segment.tau_cdag);

    // ------------  Trace ratio  -------------
        // FIXME : here we will need the K function integral 
    double ln_trace_ratio = 0;
    for (int c : seglists) {
      // Subtract the overlap with the antisegment we are removing 
      if (c != color) ln_trace_ratio = overlap(segment_t{right_segment.tau_c,left_segment.tau_cdag}, config.seglists[c]);
    }
    double trace_ratio = std::exp(ln_trace_ratio);

    // ------------  Det ratio  ---------------

    // FIXME

    // ------------  Proposition ratio ------------

    double future_number_segments = sl.size() - 1;
    double current_number_pairs = (sl.size() == 2) ? 1 : sl.size(); 
    qmc_time_t l = left_segment.tau_c - right_segment.tau_cdag; // Length of future segment 
    double prop_ratio = (future_number_segments * l * l) / current_number_pairs;

    SPDLOG_LOGGER_TRACE("trace_ratio  = {}, prop_ratio = {}, det_ratio = {}", trace_ratio, prop_ratio, det_ratio);

    return trace_ratio * det_ratio * prop_ratio;
  }

  //--------------------------------------------------

  double regroup_segment::accept() {

    SPDLOG_LOGGER_TRACE("\n - - - - - ====> ACCEPT - - - - - - - - - - -\n");

    data.dets[color].complete_operation();
    // Regroup segments 
    auto &sl = config.seglists[color];
    segment_t new_segment = segment_t{proposed_segment.tau_c,right_segment.tau_cdag}; 
    sl[left_segment_index] = new_segment; 
    sl.erase(right_segment_index);
    

    // FIXME ??? SIGNE ???
    double sign_ratio = 1; 
    return sign_ratio;
  }

  //--------------------------------------------------
  void regroup_segment::reject() {
    SPDLOG_LOGGER_TRACE("\n - - - - - ====> REJECT - - - - - - - - - - -\n");
    data.dets[color].reject_last_try();
  }
};
} // namespace moves
