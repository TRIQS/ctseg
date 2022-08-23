#include "./sperp_tau2.hpp"
#include "../logs.hpp"

namespace measures {

  sperp_tau2::sperp_tau2(params_t const &p, work_data_t const &wdata, configuration_t const &config, results_t &results)
     : wdata{wdata}, config{config}, results{results} {

    beta = p.beta;

    n_color = config.n_color();

    ss_tau2   = gf<imtime>({beta, Boson, p.n_tau_k}, {1, 1});
    ss_tau2_rev   = gf<imtime>({beta, Boson, p.n_tau_k}, {1, 1});
    ss_tau2() = 0;
    ss_tau2_rev() = 0;
    Z = 0;
  }

  // -------------------------------------

  void sperp_tau2::accumulate(double s) {

    LOG("\n =================== MEASURE < S_x S_x > (tau) = 0.5*<S_plus S_minus + S_minus S_plus> (tau) with N^2 terms ================ \n");

    double ratio = 0.0;
    Z += s; // first term in denominator coming from original segments

    for (auto const &[k1, line1] : itertools::enumerate(config.Jperp_list)) {
      
      // compute terms in numerator coming from original segments
      auto dtau11 = double(line1.tau_Splus - line1.tau_Sminus);
      auto r11 = real(wdata.Jperp(dtau11)(0, 0));
      ss_tau2[closest_mesh_pt(dtau11)] +=  0.5*s/r11;

      for (auto const &[k2, line2] : itertools::enumerate(config.Jperp_list)) {
            
        // now compute terms coming from reimagining the segments     
        if(k1 < k2) {

          // first compute the reimagined terms in the denominator
          auto dtau22 = double(line2.tau_Splus - line2.tau_Sminus);
          auto dtau12 = double(line1.tau_Splus - line2.tau_Sminus);
          auto dtau21 = double(line2.tau_Splus - line1.tau_Sminus);
          auto r22 = real(wdata.Jperp(dtau22)(0, 0));
          auto r12 = real(wdata.Jperp(dtau12)(0, 0));
          auto r21 = real(wdata.Jperp(dtau21)(0, 0));
          ratio = r12*r21/(r11*r22);
          Z += s*ratio;

          // now compute the reimagined terms in the numerator
          ss_tau2[closest_mesh_pt(dtau12)] += 0.5*s*ratio/r12;
          ss_tau2[closest_mesh_pt(dtau21)] += 0.5*s*ratio/r21;
        
          for (auto const &[k3, line3] : itertools::enumerate(config.Jperp_list)) {

            if(k3 != k1 && k3 != k2) {

              auto dtau33 = double(line3.tau_Splus - line3.tau_Sminus);
              ss_tau2[closest_mesh_pt(dtau33)] +=  0.5*s*ratio/(real(wdata.Jperp(dtau33)(0, 0)));

            }
          }
        }
      }
    }
  }   

  // -------------------------------------

  void sperp_tau2::collect_results(mpi::communicator const &c) {

    Z = mpi::all_reduce(Z, c);

    ss_tau2 = mpi::all_reduce(ss_tau2, c);

    // symmetrize and normalize
    nda::clef::placeholder<0> tau_;
    ss_tau2_rev(tau_) << ss_tau2(beta - tau_);
    ss_tau2 = 0.5*(ss_tau2 + ss_tau2_rev);
    ss_tau2 = ss_tau2 / (-beta * Z * ss_tau2.mesh().delta());

    // Fix the point at zero and beta
    ss_tau2[0] *= 2;
    ss_tau2[ss_tau2.mesh().size() - 1] *= 2;

    // store the result (not reused later, hence we can move it).
    results.sperp_tau2 = std::move(ss_tau2);
  }

} // namespace measures
