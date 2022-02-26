/*******************************************************************************
 * CTSeg: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2013-2018 by T. Ayral, H. Hafermann, P. Delange, M. Ferrero, O.
 *Parcollet Copyright (C) 2019 Simons Foundation author: N. Wentzell
 *
 * CTSEG is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * CTSEG is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * CTSEG. If not, see <http://www.gnu.org/licenses/>.
 ******************************************************************************/
#pragma once

#include <optional>
#include "types.hpp"
#include "params.hpp"
#include "work_data.hpp"
#include "inputs.hpp"
#include "results.hpp"

//#include "block_matrix.hpp"
//#include <triqs/mesh.hpp>


/// Main solver class
class solver_core {

  double beta;

  // Keep the construction parameter
  constr_params_t constr_params;

  // Keep the last solver parameters from the last call
  std::optional<solve_params_t> last_solve_params;

  //
  inputs_t inputs;
  
  // The set of results. Will be passed to measure and init by them
  results_t results;

  mpi::communicator c;
  // FIXME  Why here ??
  //double percent_done_, average_sign_;

  public:
  CPP2PY_ARG_AS_DICT solver_core(constr_params_t const &p);

  CPP2PY_ARG_AS_DICT void solve(solve_params_t const &p);

  // Green's function views for Python interface
  // do NOT add const here : python use a non const object and non const view
  // FIXME : probably fixed now, no ?

  // FIXME : Python. ?

  /// Hybridization function $\Delta^\sigma_{ab}(\tau)$
  block_gf_view<imtime> Delta_tau() { return inputs.delta; }

  /// Dynamical spin-spin interactions $\mathcal{J}_\perp(\tau)$
  gf_view<imtime> Jperp_tau() { return inputs.jperpt; }

  /// Dynamical spin-spin interactions $\mathcal{J}_\perp(\tau)$
  gf_view<imtime> D0_tau() { return inputs.d0t; }

  /// Monte Carlo sign
  //double average_sign() { return average_sign_; }
};
