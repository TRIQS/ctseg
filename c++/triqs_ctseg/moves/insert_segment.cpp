#include "insert_segment.hpp"
namespace moves {

  // FIXME : LIGNE PLEINE ?? dans overlap ???

  double insert_segment::attempt() {

    SPDLOG_LOGGER_TRACE("\n =================== ATTEMPT INSERT ================ \n");

    // --------  Selection of the times and color for the segment insertion

    // First case : there are segments. We pick up one segment, and insert at its right, i.e. smaller time

    // pick up color and first time
    color = rng(data.n_color);
    SPDLOG_LOGGER_TRACE("Inserting at color {}", n_color);

    // -------- Search for insertion window [tau1. tau2]

    qmc_time_t tau1, tau2 = data.tau_seg.get_upper_pt();

    if (not config.seglist[color].empty()) {
      // choose one segment randomly
      auto &sl       = config.seglist[color];
      long n_segment = rng(sl.size());
      // tau1 : position of cdag
      // tau2 : position of next c, possibly cyclic
      tau1                 = sl[n_segment].tau_cdag;
      bool is_last_segment = n_segment == sl.size() - 1;
      tau2                 = sl[is_last_segment ? 0 : n_segment + 1].tau_c;
    }

    // -------- Choose the new segment
    //
    auto l = tau1 - tau2;
    auto dt1     = time_point_factory.get_random_pt(l);
    auto dt2     = time_point_factory.get_random_pt(l);
    if (dt1 == dt2) return 0;
    if (dt1 < dt2) std::swap(dt1, dt2);

    proposed_segment = segment_t{tau1 + dt1, tau1 + dt2};

    SPDLOG_LOGGER_TRACE("Inserting c at{}, cdag at {}", proposed_segment.tau_c, proposed_segment.tau_cdag);

    // ----    compute ratio of trace and det  ------------

    double ln_trace_ratio = 0;
    for (auto c : range(seglists.size()))
      if (c != color) ln_trace_ratio += overlap(proposed_segment, config.seglists[c]);

    double trace_ratio = std::exp(ln_trace_ratio);

    // ----    compute ratio det  ------------

    // Try to insert in the determinants
    // FIXME : color,  seg.tau_c, seg.tau_cdag.
    double det_ratio = config->hyb_dets.try_add(seg.l, seg.r);

    // ----    compute proposition ratio ------------

    double current_number_segments = config.seglist[color];
    double prop_ratio              = (current_number_segments + 1) / (current_number_segments * l * l);

    SPDLOG_LOGGER_TRACE("trace_ratio  = {}, prop_ratio = {}, det_ratio = {}", trace_ratio, prop_ratio, det_ratio);

    return trace_ratio * det_ratio * prop_ratio;
  }

  //--------------------------------------------------

  double insert_segment::accept() {

    SPDLOG_LOGGER_TRACE("\n - - - - - ====> ACCEPT - - - - - - - - - - -\n");

    data.dets[color].complete();

    // insert the segment in an ordered list
    auto &V   = config.seglists[color];
    auto comp = [](segment_t const &s1, segment_t const &s2) { return s1.tau_c > s2.tau_c; };
     
    V.insert(pos, proposed_segment);

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
    data.dets[color].reject();

    // SPDLOG_LOGGER_TRACE("Configuration {}", config);
    // Check invariant ??
    // config->trace.check_overlap_matrix_from_scratch();
  }
};
} // namespace moves
