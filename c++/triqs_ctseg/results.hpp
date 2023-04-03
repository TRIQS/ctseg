#pragma once
#include <optional>
#include <triqs/stat/histograms.hpp>
#include <triqs/gfs.hpp>

using namespace triqs::gfs;

// Gather all the results on the CTQMC
struct results_t {

  /// Single-particle Green's function :math:`G(\tau)` in imaginary time.
  //std::optional<block_gf<imtime>> g_tau;
  block_gf<imtime> G_tau;

  /// Dynamical interaction kernel K(tau)
  gf<imtime> K_tau;

  /// Dynamical interaction kernel Kprime(tau)
  gf<imtime> Kprime_tau;

  /// Single-particle Green's function :math:`F(\tau)` in imaginary time.
  std::optional<block_gf<imtime>> F_tau;

  /// <n_a(tau) n_b(0)>
  std::optional<gf<imtime>> nn_tau;

  /// <s_x(tau) s_x(0)>
  std::optional<gf<imtime>> sperp_tau;

  /// <n_a n_b>
  std::optional<nda::matrix<double>> nn_static;

  /// Density per color. FIXME : optional ??
  nda::array<double, 1> densities;

  /// Perturbation order histogram
  std::optional<triqs::stat::histogram> perturbation_order_histo_Delta;

  /// Perturbation order histogram
  std::optional<triqs::stat::histogram> perturbation_order_histo_Jperp;

  /// Average sign 
  double sign;
};

/// writes all containers to hdf5 file
void h5_write(h5::group h5group, std::string subgroup_name, results_t const &c);

/// reads all containers to hdf5 file
void h5_read(h5::group h5group, std::string subgroup_name, results_t &c);
