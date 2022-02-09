#include "insert_segment.hpp"
namespace moves {

  // FIXME : LIGNE PLEINE ?? dans overlap ???

  double insert_segment::attempt() {

    SPDLOG_LOGGER_TRACE("\n =================== ATTEMPT INSERT ================ \n");

    // ------------ Choice of segment --------------

    // Select insertion color
    color = rng(data.n_color);
    auto &sl       = config.seglist[color];
    SPDLOG_LOGGER_TRACE("Inserting at color {}", n_color);

    // Select insertion window [tau1,tau2] (defaults to [beta,0])
    qmc_time_t tau1 = time_point_factory.get_upper_pt();
    qmc_time_t tau2 = time_point_factory.get_lower_pt();
    bool config_is_empty = sl.empty()
    if (not config_is_empty) {
      // Randomly choose one existing segment
      long ind_segment = rng(sl.size());
      tau1                 = sl[ind_segment].tau_cdag; // tau1 : pos of cdag
      bool is_last_segment = ind_segment == sl.size() - 1;
      tau2                 = sl[is_last_segment ? 0 : ind_segment + 1].tau_c; // tau2 : pos of next c, possibly cyclic
    }

    // Choose new segment within insertion window
    auto l   = tau1 - tau2;
    auto dt1 = time_point_factory.get_random_pt(l);
    auto dt2 = time_point_factory.get_random_pt(l);
    if (dt1 == dt2) return 0;
    if (dt1 > dt2) std::swap(dt1, dt2);
    proposed_segment = segment_t{tau1 - dt1, tau1 - dt2};
    // The index of the segment if it is inserted in the list of segments.
    proposed_segment_insert_pos = std::upper_bound(sl.begin(), sl.end(), proposed_segment);

    SPDLOG_LOGGER_TRACE("Inserting c at{}, cdag at {}", proposed_segment.tau_c, proposed_segment.tau_cdag);

    // ------------  Trace ratio  -------------
        // FIXME : here we will need the K function integral 
    double ln_trace_ratio = 0;
    for (int c : seglists) {
      if (c != color) ln_trace_ratio += overlap(proposed_segment, config.seglists[c]);
    }
    double trace_ratio = std::exp(ln_trace_ratio);

    // ------------  Det ratio  ---------------

    // Try insert in the det. We put all the segment in order in the dets.
    // Returns the ratio of dets (Cf det_manip doc).
    // pos is the index of the new segment, if inserted in the list.
    // FIXME : std::distance of an empoty?
    //long pos  = (config_is_empty ? 0 : std::distance(proposed_segment_insert_pos, V.begin()));
    long pos  = std::distance(proposed_segment_insert_pos, sl.begin());
    auto det_ratio = dets[color].try_insert(pos, pos, {proposed_segment.tau_c, 0}, {proposed_segment.tau_cdag, 0});

    // ------------  Proposition ratio ------------

    double current_number_segments = std::max(sl.size(),1);
    double future_number_segments = sl.size() + 1; 
    double prop_ratio              = future_number_segments / (current_number_segments * l * l);

    SPDLOG_LOGGER_TRACE("trace_ratio  = {}, prop_ratio = {}, det_ratio = {}", trace_ratio, prop_ratio, det_ratio);

    return trace_ratio * det_ratio * prop_ratio;
  }

  //--------------------------------------------------

  double insert_segment::accept() {

    SPDLOG_LOGGER_TRACE("\n - - - - - ====> ACCEPT - - - - - - - - - - -\n");

    data.dets[color].complete_operation();
    // Insert the segment in an ordered list
    config.seglists[color].insert(proposed_segment_insert_pos, proposed_segment);

    // FIXME ??? SIGNE ???
    double sign_ratio = 1; // ???config->trace.complete_insert_segment();

    // SPDLOG_LOGGER_TRACE("Configuration {}", config);
    // Check invariant ??
    // config->trace.check_overlap_matrix_from_scratch();

    return sign_ratio;
  }

  //--------------------------------------------------
  void insert_segment::reject() {
    SPDLOG_LOGGER_TRACE("\n - - - - - ====> REJECT - - - - - - - - - - -\n");
    data.dets[color].reject_last_try();

    // SPDLOG_LOGGER_TRACE("Configuration {}", config);
    // Check invariant ??
    // config->trace.check_overlap_matrix_from_scratch();
  }
};
} // namespace moves
