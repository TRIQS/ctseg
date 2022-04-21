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
    prop_ratio     = 1;

    // ----------- Propose move in each color ----------

    std::tie(idx_c_up, idx_cdag_down, tau_up)   = propose(0); // spin up
    std::tie(idx_c_down, idx_cdag_up, tau_down) = propose(1); // spin down

    // ----------- Trace ratio ----------
    // Correct for the overlap between the two modified segments
    auto &sl_up       = config.seglists[0];
    auto &sl_down     = config.seglists[1];
    auto old_seg_up   = sl_up[idx_c_up];
    auto old_seg_down = sl_down[idx_c_down];
    auto new_seg_up   = segment_t{tau_up, old_seg_up.tau_cdag};
    auto new_seg_down = segment_t{tau_down, old_seg_down.tau_cdag};
    ln_trace_ratio += -wdata.U(0, 1)
       * (overlap(new_seg_up, new_seg_down) + overlap(old_seg_up, old_seg_down) - overlap(new_seg_up, old_seg_down)
          - overlap(new_seg_down, old_seg_up));

    // Correct for the dynamical interaction between the two operators that have been moved
    if (wdata.has_Dt) {
      ln_trace_ratio -= real(wdata.K(double(tau_up - sl_down[idx_c_down].tau_c))(0, 1));
      ln_trace_ratio -= real(wdata.K(double(tau_down - sl_up[idx_c_up].tau_c))(0, 1));
      ln_trace_ratio += real(wdata.K(double(tau_down - tau_up))(0, 1));
      ln_trace_ratio += real(wdata.K(double(sl_up[idx_c_up].tau_c - sl_down[idx_c_down].tau_c))(0, 1));
    }
    double trace_ratio = std::exp(ln_trace_ratio);
    trace_ratio /= -real(wdata.Jperp(double(line.tau_Splus - line.tau_Sminus))(0, 0)) / 2;
    prop_ratio *= jl.size();

    // ----------- Det ratio -----------
    double det_ratio = 1;

    // Spin up
    auto &D_up = wdata.dets[0];
    det_ratio *= D_up.try_insert(det_lower_bound_x(D_up, sl_up[idx_cdag_up].tau_cdag), //
                                 det_lower_bound_y(D_up, tau_up),                      //
                                 {sl_up[idx_cdag_up].tau_cdag, 0}, {tau_up, 0});

    // Spin down
    auto &D_down             = wdata.dets[1];
    det_ratio *= D_down.try_insert(det_lower_bound_x(D_down, sl_down[idx_cdag_down].tau_cdag),
                                   det_lower_bound_y(D_down, tau_down), //
                                   {sl_down[idx_cdag_down].tau_cdag, 0}, {tau_down, 0});

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
    auto &sl_up               = config.seglists[0];
    auto &sl_down             = config.seglists[1];
    auto old_seg_up           = sl_up[idx_c_up];
    auto old_seg_down         = sl_down[idx_c_down];
    sl_up[idx_c_up].tau_c     = tau_up;
    sl_down[idx_c_down].tau_c = tau_down;

    // Update spin tags
    sl_up[idx_c_up].J_c           = false;
    sl_up[idx_cdag_up].J_cdag     = false;
    sl_down[idx_c_down].J_c       = false;
    sl_down[idx_cdag_down].J_cdag = false;

    // Fix segment ordering
    auto new_seg_up = sl_up[idx_c_up];
    if (is_cyclic(new_seg_up) and !is_cyclic(old_seg_up)) {
      sl_up.erase(sl_up.begin() + idx_c_up);
      sl_up.insert(sl_up.end(), new_seg_up);
    }
    if (!is_cyclic(new_seg_up) and is_cyclic(old_seg_up)) {
      sl_up.erase(sl_up.begin() + idx_c_up);
      sl_up.insert(sl_up.begin(), new_seg_up);
    }
    auto new_seg_down = sl_down[idx_c_down];
    if (is_cyclic(new_seg_down) and !is_cyclic(old_seg_down)) {
      sl_down.erase(sl_down.begin() + idx_c_down);
      sl_down.insert(sl_down.end(), new_seg_down);
    }
    if (!is_cyclic(new_seg_down) and is_cyclic(old_seg_down)) {
      sl_down.erase(sl_down.begin() + idx_c_down);
      sl_down.insert(sl_down.begin(), new_seg_down);
    }

    // Remove Jperp line
    auto &jl = config.Jperp_list;
    jl.erase(jl.begin() + line_idx);

    // Compute sign
    double final_sign = config_sign(config, wdata.dets);
    double sign_ratio = final_sign / initial_sign;
    LOG("Final sign is {}", final_sign);

    // Check invariant
    if constexpr (print_logs or ctseg_debug) check_invariant(config, wdata.dets);
    ALWAYS_EXPECTS((sign_ratio * det_sign == 1.0),
                   "Error: move has produced negative sign! Det sign is {} and additional sign is {}. Config: ",
                   det_sign, sign_ratio, config);
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

    auto &line       = config.Jperp_list[line_idx];
    int other_color  = 1 - color;
    std::string spin = (color == 0) ? "up" : "down";
    auto &sl         = config.seglists[color];
    auto &dsl        = config.seglists[other_color];

    // --------- Find the c connected to the chosen spin line ----------

    // In spin up   color, the c connected to the J line is a at tau_Sminus
    // In spin down color, the c connected to the J line is a at tau_Splus
    long idx_c  = lower_bound(sl, (color == 0 ? line.tau_Sminus : line.tau_Splus)) - sl.cbegin();
    tau_t tau_c = sl[idx_c].tau_c;

    // ---------- Find the cdag in opposite color -----------
    auto idx_cdag = cdag_in_window(tau_c + tau_t::epsilon(), tau_c - tau_t::epsilon(), dsl).back();

    // -------- Propose new position for the c ---------
    // Determine window in which the c will be moved

    auto idx_left    = (idx_c == 0) ? sl.size() - 1 : idx_c - 1;
    tau_t wtau_left  = sl[idx_left].tau_cdag;
    tau_t wtau_right = sl[idx_c].tau_cdag;
    if (idx_c == idx_left) {
      wtau_left  = tau_t::beta();
      wtau_right = tau_t::zero();
    }
    tau_t window_length = wtau_left - wtau_right;
    // Choose random time in window
    tau_t dt        = tau_t::random(rng, window_length);
    tau_t tau_c_new = sl[idx_left].tau_cdag - dt;

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
    prop_ratio *= window_length / (double(sl.size()) * cdag_in_window(wtau_left, wtau_right, dsl).size());

    return std::make_tuple(idx_c, idx_cdag, tau_c_new);
  }
}; // namespace moves
