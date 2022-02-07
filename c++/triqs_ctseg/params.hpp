#pragma once

#include "./types.hpp"

namespace triqs_ctseg {

  struct constr_params_t {

    /// Inverse temperature
    double beta;
    /// Structure of the GF (names, sizes of blocks)
    gf_struct_t gf_struct;

    /// Number of time slices for $Delta(\tau)$/$G(\tau)$/$F(\tau)$
    int n_tau = 10001;
    /// Number of time slices for $K(\tau)$
    int n_tau_k = 10001;
    /// Number of time slices for $J_\perp(\tau)$
    int n_tau_jperp = 10001;
    /// Number of Legendre coefficients for G(l)
    int n_tau_nn = 101;
    /// Number of bosonic Matsubara frequencies for $nn(i\omega)$,
    /// $\mathcal{D}_0(i\omega)$, $J_\perp(i\omega)$
    int n_w_b_nn = 32;
    /// Number of fermionic Matsubara frequencies for $G_0(i\omega)$, $G$, $F$,
    /// $\Sigma$
    int n_iw         = 1025;
    int n_legendre_g = 256;
    /// Number of time slices for $nn(\tau)$

    constr_params_t(){};
  };

  struct solve_params_t {

    /// local Hamiltonian
    Op h_int;

    /// Number of QMC cycles
    int n_cycles;

    /// Length of a single QMC cycle
    int length_cycle = 50;
    /// Number of cycles for thermalization
    int n_warmup_cycles = 5000;
    /// Seed for random number generator
    int random_seed = 34788 + 928374 * mpi::communicator().rank();
    /// Name of random number generator
    std::string random_name = "";
    /// Maximum runtime in seconds, use -1 to set infinite
    int max_time = -1;
    /// Verbosity level
    int verbosity = mpi::communicator().rank() == 0 ? 3 : 0;
    /// Whether to perform the move insert segment (see [[move_insert_segment]])
    bool move_insert_segment = true;
    /// Whether to perform the move remove segment (see [[move_remove_segment]])
    bool move_remove_segment = true;
    /// Whether to perform the move move segment (see [[move_move]])
    bool move_move = false;
    /// Whether to perform the move swap empty lines (see
    /// [[move_swap_empty_lines]])
    bool move_swap_empty_lines = false;
    /// Whether to perform the move group into spin segment (see
    /// [[move_group_into_spin_segment]])
    bool move_group_into_spin_segment = false;
    /// Whether to perform the move split spin segment (see
    /// [[move_split_spin_segment]])
    bool move_split_spin_segment = false;
    /// Whether to perform the move group into spin segment (see
    /// [[move_group_into_spin_segment2]])
    bool move_group_into_spin_segment2 = false;
    /// Whether to perform the move split spin segment (see
    /// [[move_split_spin_segment2]])
    bool move_split_spin_segment2 = false;
    /// Whether to keep Jperp negative
    bool keep_Jperp_negative = true;
    /// Whether to measure G(tau) (see [[measure_gt]])
    bool measure_gt = true;
    /// Whether to measure MC sign (see [[measure_sign]])
    bool measure_sign = true;
    /// Whether to measure improved estimator F(tau) (see [[measure_gt]]) (only
    /// isotropic spin-spin interactions if any)
    bool measure_ft = false;
    /// Whether to measure G(l) (Legendre) (see [[measure_gl]])
    bool measure_gl = false;
    /// Whether to measure improved estimator F(l) (Legendre) (see [[measure_gl]])
    /// (only isotropic spin-spin interactions if any)
    bool measure_fl = false;
    /// Whether to measure G(iomega) (see [[measure_gw]])
    bool measure_gw = false;
    /// Whether to use NFFT in the measurement of gw
    bool use_nfft_for_gw = false;
    /// Whether to measure improved estimator F(iomega) (see [[measure_gw]])(only
    /// isotropic spin-spin interactions if any)
    bool measure_fw = false;
    /// Whether to measure two-frequency correlation function (see
    /// [[measure_g2w]])
    bool measure_g2w = false;
    /// Whether to use NFFT for the precomputation of Mw
    bool use_nfft_for_Mw = false;
    /// Whether to measure two-frequency improved estimator (see [[measure_g2w]])
    bool measure_f2w = false;
    /// Whether to measure three-frequency correlation function (see
    /// [[measure_g3w]])
    bool measure_g3w = false;
    /// Whether to measure three-frequency improved estimator (see
    /// [[measure_g3w]])
    bool measure_f3w = false;
    /// Whether to measure chi_{+-}(tau) (see [[measure_chipmt]])
    bool measure_chipmt = false;
    /// Whether to measure <nn> (see [[measure_nn]])
    bool measure_nn = false;
    /// Whether to measure langle n(tau)n(0)rangle (see [[measure_nnt]])
    bool measure_nnt = false;
    /// Whether to measure chi(iomega) (see [[measure_nnw]])
    bool measure_nnw = false;
    /// Whether to evaluate vertex functions (see [[evaluate_3w_vertex]])
    bool evaluate_vertex = false;
    /// Whether to measure the perturbation order histogram (see [[measure_hist]])
    bool measure_hist = false;
    /// Whether to measure the perturbation order histogram for bosonic lines (see
    /// [[measure_hist_composite]])
    bool measure_hist_composite = false;
    /// Whether to measure histogram of impurity states (see
    /// [[measure_statehist]])
    bool measure_statehist = false;
    /// Number of fermionic Matsubara frequencies for vertex functions
    int n_w_f_vertex = 10;
    /// Number of bosonic Matsubara frequencies for vertex functions
    int n_w_b_vertex = 10;
    /// File name for 4-leg vertex gammaw
    std::string fname_gammaw = "gammaw.h5";
#ifdef HAVE_NFFT
    /// the direct algorithm is generally faster for small problem sizes, see the
    /// NFFT documentation Minimal perturbation order from which to use NFFT
    int nfft_threshold = 7;
#else /// not to have scripts fail if code is compiled without NFFT support
    /// Warning: parameter nfft_threshold will not be used because the code has
    /// been compiled without NFFT support
    int nfft_threshold = 0;
#endif
    /// Hartree shift of the chem pot
    std::vector<double> hartree_shift = std::vector<double>{};

    solve_params_t(){};
  };

  /// A struct combining both constr_params_t and solve_params_t
  struct params_t : constr_params_t, solve_params_t {
    params_t(constr_params_t const &constr_params_, solve_params_t const &solve_params_)
       : constr_params_t(constr_params_), solve_params_t(solve_params_) {}
  };

} // namespace triqs_ctseg
