#pragma once
#include <optional>

// One-particle Green's function types
using G_tau_t          = block_gf<imtime, matrix_valued>;
using G_iw_t           = block_gf<imfreq, matrix_valued>;

// Containers for measurements
struct result_set_t {

  /// Single-particle Green's function :math:`G(\tau)` in imaginary time.
  std::optional<G_tau_t> G_tau;

  /// Single-particle Green's function :math:`G(\tau)` in imaginary time.
  std::optional<G_tau_t> chi_tau;

  /// Density per color. FIXME : optional ?? 
  nda::array<double, 1> densities;

}; 

/// writes all containers to hdf5 file
void h5_write(h5::group h5group, std::string subgroup_name, container_set_t const &c);

/// reads all containers to hdf5 file
void h5_read(h5::group h5group, std::string subgroup_name, container_set_t &c);
