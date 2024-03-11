#include "invariants.hpp"
#include "logs.hpp"
#include "util.hpp"

void check_invariant(configuration_t const &config, work_data_t const &wdata) {
  check_segments(config);
  check_dets(config, wdata);
  check_jlines(config);
}

void check_segments(configuration_t const &config) {
  for (auto const &[c, sl] : itertools::enumerate(config.seglists)) {
    if (not sl.empty()) {
      for (int i = 0; i < sl.size(); ++i) {
        if (i != sl.size() - 1) {
          ALWAYS_EXPECTS((sl[i].tau_cdag > sl[i + 1].tau_c),
                         "Time order error in color {} at position {} in config \n{}", c, i, config);
          ALWAYS_EXPECTS(not is_cyclic(sl[i]), "Segment in color {} at position {} should not by cyclic in config \n{}",
                         c, i, config); // only last segment can be cyclic
        }
        ALWAYS_EXPECTS((sl[i].tau_cdag != sl[i].tau_c), "Degenerate segment in color {} at position {} in config \n{}",
                       c, i, config);
      }
    }
  }
  LOG("Segments OK.");
}

void check_dets(configuration_t const &config, work_data_t const &wdata) {
  for (auto bl : range(wdata.dets.size())) {
    auto const &D  = wdata.dets[bl];
    auto const n_orb = wdata.gf_struct[bl].second;
    // Times in det must be ordered
    if (D.size() != 0) {
      for (int i = 0; i < D.size() - 1; ++i) {
        ALWAYS_EXPECTS((D.get_x(i).first < D.get_x(i + 1).first),
                       "Det time order error (cdag) in block {} at position {} in config \n{}", bl, i, config);
        ALWAYS_EXPECTS((D.get_y(i).first < D.get_y(i + 1).first),
                       "Det time order error (c) in block {} at position {} in config \n{}", bl, i, config);
      }
    }
    // Each time in det must correspond to a time in a segment
    long n_hyb_c = 0, n_hyb_cdag = 0;
    for (auto c : range(n_orb)) {
      auto const &sl = config.seglists[wdata.block_to_color(bl, c)];
      if (not sl.empty()) {
        long det_index_c = 0, det_index_cdag = 0;
        for (auto const &seg : sl) {
          if (not seg.J_c and !is_full_line(seg)) {
            ALWAYS_EXPECTS(D.size() != 0, "Det error, block {}: there is a hybridized c but det is empty. Config: {}", bl,
                          config);
            auto det_c_time = [&](long i) { return D.get_y(i).first; };
            det_index_c     = lower_bound(det_c_time, D.size(), seg.tau_c);
            ALWAYS_EXPECTS(det_c_time(det_index_c) == seg.tau_c,
                          "Det error, block {}: tau_c = {} is not in det! Config: {}", bl, seg.tau_c, config);
            ++n_hyb_c;
          }
          if (not seg.J_cdag and !is_full_line(seg)) {
            ALWAYS_EXPECTS(D.size() != 0, "Det error, block {}: there is a hybridized cdag but det is empty. Config: {}",
                          bl, config);
            auto det_cdag_time = [&](long i) { return D.get_x(i).first; };
            det_index_cdag     = lower_bound(det_cdag_time, D.size(), seg.tau_cdag);
            ALWAYS_EXPECTS(det_cdag_time(det_index_cdag) == seg.tau_cdag,
                          "Det error, block {}: tau_cdag = {} is not in det! Config: {}", bl, seg.tau_cdag, config);
            ++n_hyb_cdag;
          }
        }
      }
    }
    ALWAYS_EXPECTS(n_hyb_c == D.size(), "Det error, block {}: missing {} c times in det. Config: {}", bl,
                   D.size() - n_hyb_c, config);
    ALWAYS_EXPECTS(n_hyb_cdag == D.size(), "Det error, block {}: missing {} cdag times in det. Config: {}", bl,
                   D.size() - n_hyb_cdag, config);
  }
  LOG("Dets OK.");
}

void check_jlines(configuration_t const &config) {
  // Prepare lists of times in spin lines
  auto const &jl = config.Jperp_list;
  if (jl.empty()) return;
  std::vector<tau_t> Splus, Sminus;
  for (auto const &[i, line] : itertools::enumerate(config.Jperp_list)) {
    Splus.push_back(line.tau_Splus);
    Sminus.push_back(line.tau_Sminus);
  }
  std::sort(Splus.begin(), Splus.end());
  std::sort(Sminus.begin(), Sminus.end());

  for (auto const &[c, sl] : itertools::enumerate(config.seglists)) {
    auto Splus_  = Splus;
    auto Sminus_ = Sminus; // will be emptied throughout the checks
    // Spin lines: each tag has to correpond to the time in a line
    for (int i = 0; i < sl.size(); ++i) {
      if (sl[i].J_c) {
        auto splus_it  = std::lower_bound(Splus_.begin(), Splus_.end(), sl[i].tau_c);
        auto sminus_it = std::lower_bound(Sminus_.begin(), Sminus_.end(), sl[i].tau_c);
        tau_t tau_plus = tau_t::zero(), tau_minus = tau_t::zero();
        if (splus_it != Splus_.end()) tau_plus = *splus_it;
        if (sminus_it != Sminus_.end()) tau_minus = *sminus_it;
        bool time_in_jlist = sl[i].tau_c == tau_plus or sl[i].tau_c == tau_minus;
        ALWAYS_EXPECTS(
           time_in_jlist,
           "Error: the c in segment at position {} in color {} has a J flag but is not in J list. Config : \n{}", i, c,
           config);
        if (sl[i].tau_c == tau_plus)
          Splus_.erase(splus_it);
        else
          Sminus_.erase(sminus_it);
      }
      if (sl[i].J_cdag) {
        auto splus_it  = std::lower_bound(Splus_.begin(), Splus_.end(), sl[i].tau_cdag);
        auto sminus_it = std::lower_bound(Sminus_.begin(), Sminus_.end(), sl[i].tau_cdag);
        tau_t tau_plus = tau_t::zero(), tau_minus = tau_t::zero();
        if (splus_it != Splus_.end()) tau_plus = *splus_it;
        if (sminus_it != Sminus_.end()) tau_minus = *sminus_it;
        bool time_in_jlist = sl[i].tau_cdag == tau_plus or sl[i].tau_cdag == tau_minus;
        ALWAYS_EXPECTS(
           time_in_jlist,
           "Error: the cdag in segment at position {} in color {} has a J flag but is not in J list. Config : \n{}", i,
           c, config);
        if (sl[i].tau_cdag == tau_plus)
          Splus_.erase(splus_it);
        else
          Sminus_.erase(sminus_it);
      }
    }
    ALWAYS_EXPECTS(Splus_.empty() and Sminus_.empty(),
                   "Error: some spin lines do not have corresponding tags on segments. Config: \n{}", config);
  }
  LOG("J lines OK.");
}