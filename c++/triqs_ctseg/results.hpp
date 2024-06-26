#pragma once
#include <optional>
#include <triqs/stat/histograms.hpp>
#include <triqs/gfs.hpp>

using namespace triqs::gfs;

// Gather all the results of the CTQMC
struct results_t {

  /// Single-particle Green's function :math:`G(\tau)`.
  block_gf<imtime> G_tau;

  /// Self-energy improved estimator :math:`F(\tau)`.
  std::optional<block_gf<imtime>> F_tau;

  /// Density-density time correlation function :math:`\langle n_a(\tau) n_b(0) \rangle`.
  std::optional<block2_gf<imtime>> nn_tau;

  /// Perpendicular spin-spin correlation function :math:`\langle s_x(\tau) s_x(0) \rangle`.
  std::optional<gf<imtime>> sperp_tau;

  /// Density-density static correlation function :math:`\langle n_a(0) n_b(0) \rangle`.
  std::optional<nda::matrix<double>> nn_static;

  /// Density per color.
  nda::array<double, 1> densities;

  /// Delta perturbation order histogram
  std::optional<triqs::stat::histogram> pert_order_histo_Delta;

  /// J_perp perturbation order histogram
  std::optional<triqs::stat::histogram> pert_order_histo_Jperp;

  /// State histogram
  std::optional<nda::vector<double>> state_hist;

  /// Average sign
  double sign;
};

/// writes all containers to hdf5 file
void h5_write(h5::group h5group, std::string subgroup_name, results_t const &c);

/// reads all containers from hdf5 file
void h5_read(h5::group h5group, std::string subgroup_name, results_t &c);
