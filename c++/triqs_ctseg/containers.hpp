#pragma once
#include <optional>
#include "types.hpp"

// Containers for inputs
struct inputs_t { 
  
  /// Hybridization function :math:`\Delta(\tau)` in imaginary time. 
  block_gf<imtime> delta;

  /// Spin-spin interaction :math:`J_{\perp}(\tau)` in imaginary time. 
  gf<imtime, matrix_valued> jperpt; // perp spin retarded interaction kernel

  /// Retarded interaction :math:`D(\tau)` in imaginary time. 
  gf<imtime, matrix_valued> d0t;

};

// Containers for measurements
struct results_t {

  /// Single-particle Green's function :math:`G(\tau)` in imaginary time.
  std::optional<G_tau_t> G_tau{};

  /// Density correlation function :math:`\Chi(\tau)` in imaginary time.
  std::optional<G_tau_t> chi_tau{};

  /// Density per color. FIXME : optional ?? 
  nda::array<double, 1> densities{};

}; 

/// writes all containers to hdf5 file
void h5_write(h5::group h5group, std::string subgroup_name, results_t const &c);

/// reads all containers to hdf5 file
void h5_read(h5::group h5group, std::string subgroup_name, results_t &c);
