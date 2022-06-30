#include "split_spin_segment.hpp"
#include "../logs.hpp"
#include <cmath>
#include <tuple>

namespace moves {

  double split_spin_segment::attempt() {

    LOG("\n =================== ATTEMPT SPLIT SPIN ================ \n");

    // ------------ Choose a spin line --------------

    auto &jl = config.Jperp_list;

    if (jl.empty()) {
      LOG("No bosonic lines!");
      return 0;
    }

    line_idx   = rng(jl.size());
    auto &line = jl[line_idx];
    LOG("Splitting S_plus at {}, S_minus at {}", line.tau_Splus, line.tau_Sminus);

    ln_trace_ratio = 0;
    prop_ratio     = jl.size();

    // ----------- Propose move in each color ----------

    std::tie(idx_c_up, idx_cdag_dn, tau_up) = propose(0); // spin up
    std::tie(idx_c_dn, idx_cdag_up, tau_dn) = propose(1); // spin down

    // ----------- Trace ratio ----------
    // Correct for the overlap between the two modified segments
    auto &sl_up     = config.seglists[0];
    auto &sl_dn     = config.seglists[1];
    auto old_seg_up = sl_up[idx_c_up];
    auto old_seg_dn = sl_dn[idx_c_dn];
    auto new_seg_up = segment_t{tau_up, old_seg_up.tau_cdag};
    auto new_seg_dn = segment_t{tau_dn, old_seg_dn.tau_cdag};

    ln_trace_ratio += -wdata.U(0, 1)
       * (overlap(new_seg_up, new_seg_dn) + overlap(old_seg_up, old_seg_dn) - //
          overlap(new_seg_up, old_seg_dn) - overlap(new_seg_dn, old_seg_up));

    // Correct for the dynamical interaction between the two operators that have been moved
    if (wdata.has_Dt) {
      ln_trace_ratio -= real(wdata.K(double(tau_up - old_seg_dn.tau_c))(0, 1));
      ln_trace_ratio -= real(wdata.K(double(tau_dn - old_seg_up.tau_c))(0, 1));
      ln_trace_ratio += real(wdata.K(double(tau_dn - tau_up))(0, 1));
      ln_trace_ratio += real(wdata.K(double(old_seg_up.tau_c - old_seg_dn.tau_c))(0, 1));
    }

    double trace_ratio = std::exp(ln_trace_ratio);
    trace_ratio /= -real(wdata.Jperp(double(line.tau_Splus - line.tau_Sminus))(0, 0)) / 2;

    // ----------- Det ratio -----------
    double det_ratio = 1;

    // Spin up
    auto &D_up = wdata.dets[0];
    det_ratio *= D_up.try_insert(det_lower_bound_x(D_up, sl_up[idx_cdag_up].tau_cdag), //
                                 det_lower_bound_y(D_up, tau_up),                      //
                                 {sl_up[idx_cdag_up].tau_cdag, 0}, {tau_up, 0});

    // Spin down
    auto &D_dn = wdata.dets[1];
    det_ratio *=
       D_dn.try_insert(det_lower_bound_x(D_dn, sl_dn[idx_cdag_dn].tau_cdag), det_lower_bound_y(D_dn, tau_dn), //
                       {sl_dn[idx_cdag_dn].tau_cdag, 0}, {tau_dn, 0});

    LOG("trace_ratio  = {}, prop_ratio = {}, det_ratio = {}", trace_ratio, prop_ratio, det_ratio);

    double prod = trace_ratio * det_ratio * prop_ratio;
    det_sign    = (det_ratio > 0) ? 1.0 : -1.0;
    return (std::isfinite(prod) ? prod : det_sign);
  }

  //--------------------------------------------------

  double split_spin_segment::accept() {

    LOG("\n - - - - - ====> ACCEPT - - - - - - - - - - -\n");

    double initial_sign = config_sign(config, wdata.dets);
    LOG("Initial sign is {}. Initial configuration: {}", initial_sign, config);

    // Update the dets
    wdata.dets[0].complete_operation();
    wdata.dets[1].complete_operation();

    // Update the segments
    auto &sl_up = config.seglists[0];
    auto &sl_dn = config.seglists[1];

    sl_up[idx_c_up].tau_c = tau_up;
    sl_dn[idx_c_dn].tau_c = tau_dn;

    // Update spin tags
    sl_up[idx_c_up].J_c       = false;
    sl_up[idx_cdag_up].J_cdag = false;
    sl_dn[idx_c_dn].J_c       = false;
    sl_dn[idx_cdag_dn].J_cdag = false;

    fix_ordering_first_last(sl_up);
    fix_ordering_first_last(sl_dn);

    // Remove Jperp line
    auto &jl = config.Jperp_list;
    jl.erase(begin(jl) + line_idx);

    // Compute sign
    double final_sign = config_sign(config, wdata.dets);
    double sign_ratio = final_sign / initial_sign;
    LOG("Final sign is {}", final_sign);

    // Check invariant
    if constexpr (print_logs or ctseg_debug) check_invariant(config, wdata.dets);
    // ALWAYS_EXPECTS((sign_ratio * det_sign == 1.0),
    //                "Error: move has produced negative sign! Det sign is {} and additional sign is {}. Config: ",
    //                det_sign, sign_ratio, config);
    if (sign_ratio * det_sign == -1.0) spdlog::info("WARNING: move produced negative sign!");
    LOG("Configuration is {}", config);

    return sign_ratio;
  }

  //--------------------------------------------------

  void split_spin_segment::reject() {

    LOG("\n - - - - - ====> REJECT - - - - - - - - - - -\n");
    wdata.dets[0].reject_last_try();
    wdata.dets[1].reject_last_try();
  }

  //--------------------------------------------------

  std::tuple<long, long, tau_t> split_spin_segment::propose(int color) {

    auto &line      = config.Jperp_list[line_idx];
    int other_color = 1 - color;
    auto &sl        = config.seglists[color];
    auto &dsl       = config.seglists[other_color];

    // --------- Find the c connected to the chosen spin line ----------

    // In spin up   color, the c connected to the J line is a at tau_Sminus
    // In spin down color, the c connected to the J line is a at tau_Splus
    long idx_c = lower_bound(sl, (color == 0 ? line.tau_Sminus : line.tau_Splus)) - sl.cbegin();
    auto tau_c = sl[idx_c].tau_c;

    // ---------- Find the cdag in opposite color -----------

    // FIXME : ok, the vector is always of size 1 ...
    auto idx_cdag = cdag_in_window(tau_c + tau_t::epsilon(), tau_c - tau_t::epsilon(), dsl).back();

    // -------- Propose new position for the c ---------

    // Determine window in which the c will be moved
    auto idx_left   = (idx_c == 0) ? sl.size() - 1 : idx_c - 1;
    auto wtau_left  = sl[idx_left].tau_cdag;
    auto wtau_right = sl[idx_c].tau_cdag;
    if (idx_c == idx_left) { // only one segment ...
      wtau_left  = tau_t::beta();
      wtau_right = tau_t::zero();
    }
    tau_t window_length = wtau_left - wtau_right;

    // Choose random time in window left-right (possibly cyclic)
    auto dt        = tau_t::random(rng, window_length);
    auto tau_c_new = sl[idx_left].tau_cdag - dt;

    LOG("Spin {}: moving c at position {} from {} to {}.", (color == 0) ? "up" : "down", idx_c, tau_c, tau_c_new);
    auto new_seg = segment_t{tau_c_new, sl[idx_c].tau_cdag};

    // -------- Trace ratio ---------
    ln_trace_ratio += wdata.mu(color) * (double(new_seg.length()) - double(sl[idx_c].length()));
    for (auto const &[c, slc] : itertools::enumerate(config.seglists)) {
      if (c != color) {
        ln_trace_ratio += -wdata.U(c, color) * overlap(slc, new_seg);
        ln_trace_ratio -= -wdata.U(c, color) * overlap(slc, sl[idx_c]);
      }
      if (wdata.has_Dt) {
        ln_trace_ratio += K_overlap(slc, tau_c_new, true, wdata.K, c, color);
        ln_trace_ratio -= K_overlap(slc, tau_c, true, wdata.K, c, color);
      }
    }
    if (wdata.has_Dt) ln_trace_ratio -= real(wdata.K(double(tau_c_new - tau_c))(color, color));

    // --------- Prop ratio ---------
    // T direct  = 1/window_length
    // T inverse =
    prop_ratio *= window_length / (double(sl.size()) * cdag_in_window(wtau_left, wtau_right, dsl).size());

    return {idx_c, idx_cdag, tau_c_new};
  }
}; // namespace moves
