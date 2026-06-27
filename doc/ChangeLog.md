# Changelog

## Version 4.0.0

ctseg version 4.0.0 is a feature release that is compatible with TRIQS version 4.0.0
and includes an update to the latest app4triqs skeleton. The headline additions are a
new functional `solve_generic()` interface, three- and four-point vertex measurements,
a `<n(nu) n(-nu)>` (`nn_nu_dlr`) measurement, and double insert/remove segment moves
(now enabled by default). The Python bindings have been migrated from cpp2py to
clair + c2py.

We thank all contributors: Thomas Ayral, Jennifer Coulter, Thomas Hahn, Alexander Hampel,
Nikita Kavokine, Hao Lu, Erik van Loon, Henri Menke, Olivier Parcollet, Dylan Simon,
Hugo U. R. Strand, Nils Wentzell

### Porting and compatibility

Following the TRIQS 4.0 release, several Python module locations changed. In particular the
`triqs.gf` package was renamed to `triqs.gfs`; update your imports accordingly (e.g.
`from triqs.gfs import *`). Run the `port_to_triqs4` script from the TRIQS distribution to
apply the mechanical renames to your own scripts.

The Python bindings are now generated with clair + c2py instead of cpp2py. The bound
`Solver` / `solver_core` interface is unchanged.

The `nn_nu` measurement has been renamed to `nn_nu_dlr`.

### General
* Add a functional `solve_generic()` interface for ctseg
* Add three- and four-point vertex measurements (g2w, g3w) and expose them as block Green's functions (#9)
* Implement the `<n(nu) n(-nu)>` (`nn_nu_dlr`) measurement
* Implement double insert/remove segment moves and enable `move_double` by default
* Implement configuration visualization
* Handle complex-flagged `h_loc0`: take the real part within `imag_threshold` (#28)
* Fix a bug when setting detmanip parameters
* Fix `try_insert2` for off-diagonal blocks
* Migrate the Python bindings from cpp2py to clair + c2py
* Fix gcc compilation by removing the specialized formatter for nda arrays
* Fix build failure with fmt 11 by using `fmt::streamed` for nda types
* Rename the `nn_nu` measurement to `nn_nu_dlr`
* Update copyright/license headers and notices

### doc
* Add a one-boson retarded interaction tutorial
* Document the three- and four-point vertex measurements
* Add a section on tests and conditions for half-filling (#12)
* Fix Sphinx warnings and errors
* Add Flatiron Institute support notice to README.md
* Add Debian package and Anaconda installation instructions

### jenkins / cmake / ghactions
* Use the latest app4triqs skeleton (k8s Jenkins setup, doxygen update)
* Bump spdlog to 1.15.1
* Synchronize the GitHub Actions build workflow with cthyb

## Version 3.3.0

This is the initial release of ctseg, a TRIQS based segment picture impurity solver
with spin-spin interactions.  It is compatible against TRIQS Version 3.3.

We thank all contributors: Nikita Kavokine, Hao Lu, Alexander Hampel, Nils Wentzell, Olivier Parcollet
