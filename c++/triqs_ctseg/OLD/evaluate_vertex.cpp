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
#include "./evaluate_vertex.hpp"

namespace triqs_ctseg {

void evaluate_3w_vertex(block_gf<imfreq> const &gw, block_gf<imfreq> const &fw,
                        gf_3w_container_t const &g3w,
                        gf_3w_container_t const &f3w, bool measure_g3w,
                        bool measure_f3w, std::string fname) {
  mpi::communicator c;

  if (!(measure_g3w || measure_f3w))
    return;

  if (c.rank() == 0) {
    nda::clef::placeholder<0> U_;
    nda::clef::placeholder<1> V_;
    nda::clef::placeholder<2> a_;
    nda::clef::placeholder<3> b_;
    nda::clef::placeholder<4> c_;
    nda::clef::placeholder<5> d_;
    nda::clef::placeholder<6> inu_;
    nda::clef::placeholder<7> inup_;
    nda::clef::placeholder<8> iom_;
    nda::clef::placeholder<9> u_;
    nda::clef::placeholder<10> v_;
    nda::clef::placeholder<11> w_;
    nda::clef::placeholder<12> x_;
    nda::clef::placeholder<13> f_;

    double beta = (double)gw[0].mesh().domain().beta;
    int n_blocks = gw.size();

    // convert g3w to block_gf2
    auto g3w_2ind = make_block2_gf(n_blocks, n_blocks, g3w[0]);
    auto f3w_2ind = make_block2_gf(n_blocks, n_blocks, f3w[0]);
    for (int U = 0; U < n_blocks; U++)
      for (int V = 0; V < n_blocks; V++) {
        g3w_2ind(U, V) = g3w[U * n_blocks + V];
        f3w_2ind(U, V) = f3w[U * n_blocks + V];
      }

    auto gammaw = g3w_2ind;
    auto g3w_disc = g3w_2ind;
    auto g3w_conn_1 = g3w_2ind;
    auto g3w_conn_3 = g3w_2ind;
    for (auto &g : gammaw)
      g() = 0.;

    auto g3w_disconnected =
        -beta *
        ((gw[U_](inu_)(b_, a_) * gw[V_](inup_)(d_, c_) * kronecker(iom_)) -
         (gw[U_](inu_)(d_, a_) * gw[U_](inu_ + iom_)(b_, c_) *
          kronecker(inu_, inup_) * kronecker(U_, V_)));

    auto g3w_connected_1 =
        g3w_2ind(U_, V_)(inu_, inup_, iom_)(a_, b_, c_, d_) - g3w_disconnected;

    // g3w_disc(U_,V_)(inu_,inup_,iom_)(a_,b_,c_,d_) <<   g3w_disconnected;
    g3w_conn_1(U_, V_)(inu_, inup_, iom_)(a_, b_, c_, d_)
        << g3w_2ind(U_, V_)(inu_, inup_, iom_)(a_, b_, c_, d_) -
               g3w_disconnected;
    /*
    auto g3w_connected_2 =
        (gw(U_,V_)(inu_)(0, 0) * f3w(U_, V_)(inu_, inup_, iom_) -
         fw(U_,V_)(inu_)(0, 0) * g3w_disconnected) / (1.0 + fw[U_](inu_)(0, 0));
    */
    int bl_size = gw[0].target_shape()[0];
    auto g3w_connected_3 =
        sum(-gw[U_](inu_)(b_, f_) *
                    f3w_2ind(U_, V_)(inu_, inup_, iom_)(f_, a_, c_, d_) -
                fw[U_](inu_)(b_, f_) *
                    g3w_2ind(U_, V_)(inu_, inup_, iom_)(a_, f_, c_, d_),
            f_ = range(0, bl_size));
    g3w_conn_3(U_, V_)(inu_, inup_, iom_)(a_, b_, c_, d_) << g3w_connected_3;

    auto inv_gw = 1 / gw;
    auto legs = (inv_gw[U_](inu_)(b_, u_) * inv_gw[U_](inu_ + iom_)(v_, a_) *
                 inv_gw[V_](inup_ + iom_)(d_, w_) * inv_gw[V_](inup_)(x_, c_));
    if (measure_g3w)
      gammaw(U_, V_)(inu_, inup_, iom_)(u_, v_, w_, x_)
          << sum(g3w_connected_1 * legs, a_ = range(0, bl_size),
                 b_ = range(0, bl_size), c_ = range(0, bl_size),
                 d_ = range(0, bl_size));
    /*
    if (measure_f3w)
      gammaw[U_, V_](inu_, inup_, iom_) << g3w_connected_2 / legs;
    */
    if (measure_g3w && measure_f3w)
      gammaw(U_, V_)(inu_, inup_, iom_)(u_, v_, w_, x_)
          << sum(g3w_connected_3 * legs, a_ = range(0, bl_size),
                 b_ = range(0, bl_size), c_ = range(0, bl_size),
                 d_ = range(0, bl_size));

    // convert back to block_gf
    auto gammaw_1ind = g3w;
    auto g3w_conn_1_1ind = g3w;
    auto g3w_conn_3_1ind = g3w;
    // auto g3w_disc_1ind = g3w;
    for (int U = 0; U < n_blocks; U++)
      for (int V = 0; V < n_blocks; V++) {
        gammaw_1ind[U * n_blocks + V] = gammaw(U, V);
        g3w_conn_1_1ind[U * n_blocks + V] = g3w_conn_1(U, V);
        g3w_conn_3_1ind[U * n_blocks + V] = g3w_conn_3(U, V);
        // g3w_disc_1ind[U*n_blocks+V] = g3w_disc[{U,V}] ;
      }

    h5::file gamma_file(fname.c_str(), 'w');
    h5_write(gamma_file, "gammaw", gammaw_1ind);
    h5_write(gamma_file, "g3w_conn_1", g3w_conn_1_1ind);
    h5_write(gamma_file, "g3w_conn_3", g3w_conn_3_1ind);
    // h5_write(gamma_file, "g3w_disc", g3w_disc_1ind);

    gamma_file.close();
  }
}
/**
void evaluate_2w_vertex(const boost::mpi::communicator &c,
                        const qmc_parameters *params, block_gf<imfreq> &gw,
                        block_gf<imfreq> &fw, gf_2w_container_t &g2w,
                        gf_2w_container_t &f2w, gf<imfreq> &nnw,
                        matrix<double> &nn_matrix
                        ) {

  if(!(params->measure_g2w || params->measure_f2w)) return;

  if (c.rank() == 0) {
    nda::clef::placeholder<0> a_;
    nda::clef::placeholder<1> b_;
    nda::clef::placeholder<2> inu_;
    nda::clef::placeholder<4> iom_;

    double beta = (double) params->beta;

    int n_w_f_vertex = std::get<0>(g2w.mesh().components()).last_index()+1;
    int n_w_b = std::get<1>(g2w.mesh().components()).last_index()+1;

    array<gf_2w_container_t, 1> lambda_charge(params->n_color);
    lambda_charge() = gf_2w_container_t {
      { { beta, Fermion, n_w_f_vertex, imfreq::option::all_frequencies } , {
beta, Boson, n_w_b, imfreq::option::positive_frequencies_only } }, {1,1}
    }
    ;

    for (auto &g : lambda_charge)
      g() = 0.;

    double n0 = nn_matrix(0, 0);
    double n1 = nn_matrix(1, 1);

    auto chi_charge = -(
        nnw(iom_)(0, 0) + nnw(iom_)(0, 1) + nnw(iom_)(1, 0) + nnw(iom_)(1, 1) -
        beta * (n0 * n0 + n0 * n1 + n1 * n0 + n1 * n1) *
            kronecker(iom_)); // later use sum(nnw, a_, b_, domain)

    auto legs = gw[a_](inu_)(0, 0) * gw[a_](inu_ + iom_)(0, 0);

    auto g2w_charge =
        eval(g2w(a_, b_)(inu_, iom_), b_ = 0) +
        eval(g2w(a_, b_)(inu_, iom_), b_ = 1); // later use sum over b_ here
    auto f2w_charge =
        eval(f2w(a_, b_)(inu_, iom_), b_ = 0) +
        eval(f2w(a_, b_)(inu_, iom_), b_ = 1); // later use sum over b_ here

    auto g2w_disconnected =
        beta * gw[a_](inu_)(0, 0) * (n0 + n1) * kronecker(iom_) -
        gw[a_](inu_)(0, 0) * gw[a_](inu_ + iom_)(0, 0); //sum over flavors of
the density in the first term

    auto g2w_connected_1 = g2w_charge - g2w_disconnected; //if only g is known
    auto g2w_connected_2 = (gw[a_](inu_)(0, 0) * f2w_charge -
                            fw[a_](inu_)(0, 0) * g2w_disconnected) /
                           (1.0 + fw[a_](inu_)(0, 0));    //if only f is known
    auto g2w_connected_3 =
        gw[a_](inu_)(0, 0) * f2w_charge -
        fw[a_](inu_)(0, 0) * g2w_charge;                  //if g and f are known

    if (params->measure_g2w) lambda_charge(a_)(inu_, iom_) << (g2w_connected_1 /
legs - 1.0) / chi_charge; if (params->measure_f2w) lambda_charge(a_)(inu_, iom_)
<< (g2w_connected_2 / legs - 1.0) / chi_charge; if (params->measure_g2w &&
params->measure_f2w) lambda_charge(a_)(inu_, iom_) << (g2w_connected_3 / legs
- 1.0) / chi_charge;

    h5::file lambda_file(params->fname_lambdaw.c_str(), 'w');
    h5_write(lambda_file, "lambdaw_charge", lambda_charge);

  }
}
*/

} // namespace triqs_ctseg
