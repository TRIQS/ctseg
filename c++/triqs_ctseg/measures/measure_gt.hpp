/*******************************************************************************
 * CTSEG: TRIQS hybridization-expansion segment solver
 *
 * Copyright (C) 2013-2018 by T. Ayral, H. Hafermann, P. Delange, M. Ferrero, O.
 *Parcollet
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
#include "../configuration.hpp"
#include "../work_data.hpp"
#include "../results.hpp"
//#include "../precompute_fprefactor.hpp"

namespace measures {

  //using namespace triqs::gfs;
  //using namespace triqs::mesh;

  struct measure_g_f_tau {

    work_data_t const &wdata;
    configuration_t const &config;
    results_t & results;
    double beta;

    block_gf<imtime> g_tau, f_tau;

    
    // The prefactor of integrals
    //std::shared_ptr<precompute_fprefactor> fprefactor;

    double Z;

    // double beta, Noverbeta, Z;

    //accumulator<double> gt_stack = {0.0, -1, -1};
    //accumulator<double> Z_stack  = {0.0, -1, -1};
    //int counter;
    //double accum;

    measure_g_f_tau(params_t const & params, work_data_t const &wdata, configuration_t const &config, results_t &results);

    void accumulate(double s);
    void collect_results(mpi::communicator const &c);
  };
} // namespace measures
