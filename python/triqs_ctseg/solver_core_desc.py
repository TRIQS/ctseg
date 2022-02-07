# Generated automatically using the command :
# c++2py ../../c++/triqs_ctseg/solver_core.hpp -p --members_read_only -N triqs_ctseg -a triqs_ctseg -m solver_core -o solver_core --moduledoc="The python module for triqs_ctseg" -C triqs -C nda_py -C triqs_ctseg --includes="../../c++" --cxxflags="-std=c++20" --target_file_only
from cpp2py.wrap_generator import *

# The module
module = module_(full_name = "solver_core", doc = r"The python module for triqs_ctseg", app_name = "triqs_ctseg")

# Imports
module.add_imports(*['triqs_ctseg.block_matrix', 'triqs.gf', 'triqs.gf.meshes', 'triqs.operators'])

# Add here all includes
module.add_include("triqs_ctseg/solver_core.hpp")

# Add here anything to add in the C++ code at the start, e.g. namespace using
module.add_preamble("""
#include <cpp2py/converters/optional.hpp>
#include <cpp2py/converters/pair.hpp>
#include <cpp2py/converters/string.hpp>
#include <cpp2py/converters/vector.hpp>
#include <nda_py/cpp2py_converters.hpp>
#include <triqs/cpp2py_converters/gf.hpp>
#include <triqs/cpp2py_converters/operators_real_complex.hpp>
#include <triqs/cpp2py_converters/real_or_complex.hpp>

using namespace triqs_ctseg;
""")


# The class solver_core
c = class_(
        py_type = "SolverCore",  # name of the python class
        c_type = "triqs_ctseg::solver_core",   # name of the C++ class
        doc = r"""Main solver class


  Worker which runs the quantum Monte-Carlo simulation.""",   # doc of the C++ class
        hdf5 = False,
)

c.add_member(c_name = "constr_params",
             c_type = "triqs_ctseg::constr_params_t",
             read_only= True,
             doc = r"""""")

c.add_member(c_name = "last_solve_params",
             c_type = "std::optional<solve_params_t>",
             read_only= True,
             doc = r"""""")

c.add_constructor("""(**triqs_ctseg::constr_params_t)""", doc = r"""constructor


                           Allocates the main observables.

+----------------+--------------------------+---------+----------------------------------------------------------------------------------------------------------+
| Parameter Name | Type                     | Default | Documentation                                                                                            |
+================+==========================+=========+==========================================================================================================+
| beta           | double                   | --      | Inverse temperature                                                                                      |
+----------------+--------------------------+---------+----------------------------------------------------------------------------------------------------------+
| gf_struct      | triqs_ctseg::gf_struct_t | --      | Structure of the GF (names, sizes of blocks)                                                             |
+----------------+--------------------------+---------+----------------------------------------------------------------------------------------------------------+
| n_tau          | int                      | 10001   | Number of time slices for $Delta(\tau)$/$G(\tau)$/$F(\tau)$                                              |
+----------------+--------------------------+---------+----------------------------------------------------------------------------------------------------------+
| n_tau_k        | int                      | 10001   | Number of time slices for $K(\tau)$                                                                      |
+----------------+--------------------------+---------+----------------------------------------------------------------------------------------------------------+
| n_tau_jperp    | int                      | 10001   | Number of time slices for $J_\perp(\tau)$                                                                |
+----------------+--------------------------+---------+----------------------------------------------------------------------------------------------------------+
| n_tau_nn       | int                      | 101     | Number of Legendre coefficients for G(l)                                                                 |
+----------------+--------------------------+---------+----------------------------------------------------------------------------------------------------------+
| n_w_b_nn       | int                      | 32      | Number of bosonic Matsubara frequencies for $nn(i\omega)$, $\mathcal{D}_0(i\omega)$, $J_\perp(i\omega)$  |
+----------------+--------------------------+---------+----------------------------------------------------------------------------------------------------------+
| n_iw           | int                      | 1025    | Number of fermionic Matsubara frequencies for $G_0(i\omega)$, $G$, $F$, $\Sigma$                         |
+----------------+--------------------------+---------+----------------------------------------------------------------------------------------------------------+
| n_legendre_g   | int                      | 256     |                                                                                                          |
+----------------+--------------------------+---------+----------------------------------------------------------------------------------------------------------+
""")

c.add_method("""void solve (**triqs_ctseg::solve_params_t)""",
             doc = r"""solve method: starts the Metropolis algorithm


 Steps:

 - extract :math:`\Delta^\sigma_{ab}(\tau)` and :math:`\mu^\sigma_a` from
 :math:`\mathcal{G}^\sigma_{ab}(i\omega)`

 - if :math:`\mathcal{D}^{\sigma\sigma'}_{0,ab}(i\Omega)\neq 0`,  extract
 :math:`K(^{\sigma\sigma'}_{ab}\tau)` and $\partial_\tau
 K^{\sigma\sigma'}_{ab}(\tau)$ from
 :math:`\mathcal{D}^{\sigma\sigma'}_{0,ab}(i\Omega)`

 - if :math:`\mathcal{J}_{\perp,a}(i\Omega)\neq 0`,  extract $\partial_\tau
 K_{\perp,a}(\tau)$ from :math:`\mathcal{J}_{\perp,a}(i\Omega)`
 - add the moves and measures according to the parameters supplied by
 the user

 - start the Monte-Carlo simulation

 - finalize the Monte Carlo simulation

+-------------------------------+---------------------+-----------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| Parameter Name                | Type                | Default                                 | Documentation                                                                                                                                                                 |
+===============================+=====================+=========================================+===============================================================================================================================================================================+
| h_int                         | triqs_ctseg::Op     | --                                      | local Hamiltonian                                                                                                                                                             |
+-------------------------------+---------------------+-----------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| n_cycles                      | int                 | --                                      | Number of QMC cycles                                                                                                                                                          |
+-------------------------------+---------------------+-----------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| length_cycle                  | int                 | 50                                      | Length of a single QMC cycle                                                                                                                                                  |
+-------------------------------+---------------------+-----------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| n_warmup_cycles               | int                 | 5000                                    | Number of cycles for thermalization                                                                                                                                           |
+-------------------------------+---------------------+-----------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| random_seed                   | int                 | 34788+928374*mpi::communicator().rank() | Seed for random number generator                                                                                                                                              |
+-------------------------------+---------------------+-----------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| random_name                   | std::string         | ""                                      | Name of random number generator                                                                                                                                               |
+-------------------------------+---------------------+-----------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| max_time                      | int                 | -1                                      | Maximum runtime in seconds, use -1 to set infinite                                                                                                                            |
+-------------------------------+---------------------+-----------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| verbosity                     | int                 | mpi::communicator().rank()==0?3:0       | Verbosity level                                                                                                                                                               |
+-------------------------------+---------------------+-----------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| move_insert_segment           | bool                | true                                    | Whether to perform the move insert segment (see [[move_insert_segment]])                                                                                                      |
+-------------------------------+---------------------+-----------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| move_remove_segment           | bool                | true                                    | Whether to perform the move remove segment (see [[move_remove_segment]])                                                                                                      |
+-------------------------------+---------------------+-----------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| move_move                     | bool                | false                                   | Whether to perform the move move segment (see [[move_move]])                                                                                                                  |
+-------------------------------+---------------------+-----------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| move_swap_empty_lines         | bool                | false                                   | Whether to perform the move swap empty lines (see [[move_swap_empty_lines]])                                                                                                  |
+-------------------------------+---------------------+-----------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| move_group_into_spin_segment  | bool                | false                                   | Whether to perform the move group into spin segment (see [[move_group_into_spin_segment]])                                                                                    |
+-------------------------------+---------------------+-----------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| move_split_spin_segment       | bool                | false                                   | Whether to perform the move split spin segment (see [[move_split_spin_segment]])                                                                                              |
+-------------------------------+---------------------+-----------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| move_group_into_spin_segment2 | bool                | false                                   | Whether to perform the move group into spin segment (see [[move_group_into_spin_segment2]])                                                                                   |
+-------------------------------+---------------------+-----------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| move_split_spin_segment2      | bool                | false                                   | Whether to perform the move split spin segment (see [[move_split_spin_segment2]])                                                                                             |
+-------------------------------+---------------------+-----------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| keep_Jperp_negative           | bool                | true                                    | Whether to keep Jperp negative                                                                                                                                                |
+-------------------------------+---------------------+-----------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| measure_gt                    | bool                | true                                    | Whether to measure G(tau) (see [[measure_gt]])                                                                                                                                |
+-------------------------------+---------------------+-----------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| measure_sign                  | bool                | true                                    | Whether to measure MC sign (see [[measure_sign]])                                                                                                                             |
+-------------------------------+---------------------+-----------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| measure_ft                    | bool                | false                                   | Whether to measure improved estimator F(tau) (see [[measure_gt]]) (only isotropic spin-spin interactions if any)                                                              |
+-------------------------------+---------------------+-----------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| measure_gl                    | bool                | false                                   | Whether to measure G(l) (Legendre) (see [[measure_gl]])                                                                                                                       |
+-------------------------------+---------------------+-----------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| measure_fl                    | bool                | false                                   | Whether to measure improved estimator F(l) (Legendre) (see [[measure_gl]]) (only isotropic spin-spin interactions if any)                                                     |
+-------------------------------+---------------------+-----------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| measure_gw                    | bool                | false                                   | Whether to measure G(iomega) (see [[measure_gw]])                                                                                                                             |
+-------------------------------+---------------------+-----------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| use_nfft_for_gw               | bool                | false                                   | Whether to use NFFT in the measurement of gw                                                                                                                                  |
+-------------------------------+---------------------+-----------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| measure_fw                    | bool                | false                                   | Whether to measure improved estimator F(iomega) (see [[measure_gw]])(only isotropic spin-spin interactions if any)                                                            |
+-------------------------------+---------------------+-----------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| measure_g2w                   | bool                | false                                   | Whether to measure two-frequency correlation function (see [[measure_g2w]])                                                                                                   |
+-------------------------------+---------------------+-----------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| use_nfft_for_Mw               | bool                | false                                   | Whether to use NFFT for the precomputation of Mw                                                                                                                              |
+-------------------------------+---------------------+-----------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| measure_f2w                   | bool                | false                                   | Whether to measure two-frequency improved estimator (see [[measure_g2w]])                                                                                                     |
+-------------------------------+---------------------+-----------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| measure_g3w                   | bool                | false                                   | Whether to measure three-frequency correlation function (see [[measure_g3w]])                                                                                                 |
+-------------------------------+---------------------+-----------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| measure_f3w                   | bool                | false                                   | Whether to measure three-frequency improved estimator (see [[measure_g3w]])                                                                                                   |
+-------------------------------+---------------------+-----------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| measure_chipmt                | bool                | false                                   | Whether to measure chi_{+-}(tau) (see [[measure_chipmt]])                                                                                                                     |
+-------------------------------+---------------------+-----------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| measure_nn                    | bool                | false                                   | Whether to measure <nn> (see [[measure_nn]])                                                                                                                                  |
+-------------------------------+---------------------+-----------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| measure_nnt                   | bool                | false                                   | Whether to measure langle n(tau)n(0)rangle (see [[measure_nnt]])                                                                                                              |
+-------------------------------+---------------------+-----------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| measure_nnw                   | bool                | false                                   | Whether to measure chi(iomega) (see [[measure_nnw]])                                                                                                                          |
+-------------------------------+---------------------+-----------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| evaluate_vertex               | bool                | false                                   | Whether to evaluate vertex functions (see [[evaluate_3w_vertex]])                                                                                                             |
+-------------------------------+---------------------+-----------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| measure_hist                  | bool                | false                                   | Whether to measure the perturbation order histogram (see [[measure_hist]])                                                                                                    |
+-------------------------------+---------------------+-----------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| measure_hist_composite        | bool                | false                                   | Whether to measure the perturbation order histogram for bosonic lines (see [[measure_hist_composite]])                                                                        |
+-------------------------------+---------------------+-----------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| measure_statehist             | bool                | false                                   | Whether to measure histogram of impurity states (see [[measure_statehist]])                                                                                                   |
+-------------------------------+---------------------+-----------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| n_w_f_vertex                  | int                 | 10                                      | Number of fermionic Matsubara frequencies for vertex functions                                                                                                                |
+-------------------------------+---------------------+-----------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| n_w_b_vertex                  | int                 | 10                                      | Number of bosonic Matsubara frequencies for vertex functions                                                                                                                  |
+-------------------------------+---------------------+-----------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| fname_gammaw                  | std::string         | "gammaw.h5"                             | File name for 4-leg vertex gammaw                                                                                                                                             |
+-------------------------------+---------------------+-----------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| nfft_threshold                | int                 | 0                                       | not to have scripts fail if code is compiled without NFFT support Warning: parameter nfft_threshold will not be used because the code has been compiled without NFFT support  |
+-------------------------------+---------------------+-----------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| hartree_shift                 | std::vector<double> | std::vector<double>{}                   | Hartree shift of the chem pot                                                                                                                                                 |
+-------------------------------+---------------------+-----------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
""")

c.add_method("""void sanity_check (triqs_ctseg::solve_params_t p, int n_w, int n_w_b)""",
             doc = r"""""")

c.add_property(name = "G0_iw",
               getter = cfunction("block_gf_view<triqs::mesh::imfreq> G0_iw ()"),
               doc = r"""Weiss field $\mathcal{G}^{\sigma}_{0,ab}(i\omega)$""")

c.add_property(name = "D0_iw",
               getter = cfunction("block_gf_view<triqs::mesh::imfreq> D0_iw ()"),
               doc = r"""Density-density retarded interactions
   $\mathcal{D}^{\sigma\sigma'}_{0,ab}(i\Omega)$""")

c.add_property(name = "Jperp_iw",
               getter = cfunction("gf_view<triqs::mesh::imfreq> Jperp_iw ()"),
               doc = r"""Dynamical spin-spin interaction, perpendicual components:
   $\mathcal{J}_\perp(i\Omega)$""")

c.add_property(name = "Delta_tau",
               getter = cfunction("block_gf_view<triqs::mesh::imtime> Delta_tau ()"),
               doc = r"""Hybridization function $\Delta^\sigma_{ab}(\tau)$""")

c.add_property(name = "K_tau",
               getter = cfunction("block_gf_view<triqs::mesh::imtime> K_tau ()"),
               doc = r"""Dynamical kernel $K(\tau)$""")

c.add_property(name = "Kprime_tau",
               getter = cfunction("block_gf_view<triqs::mesh::imtime> Kprime_tau ()"),
               doc = r"""Derivative of the dynamical kernel $\partial_\tau K(\tau)$""")

c.add_property(name = "Kperpprime_tau",
               getter = cfunction("gf_view<triqs::mesh::imtime> Kperpprime_tau ()"),
               doc = r"""Derivative of the dynamical kernel $\partial_\tau K_\perp(\tau)$""")

c.add_property(name = "Jperp_tau",
               getter = cfunction("gf_view<triqs::mesh::imtime> Jperp_tau ()"),
               doc = r"""Dynamical spin-spin interactions $\mathcal{J}_\perp(\tau)$""")

c.add_property(name = "G_tau",
               getter = cfunction("block_gf_view<triqs::mesh::imtime> G_tau ()"),
               doc = r"""Impurity Green's function $G^\sigma_{ab}(\tau)$ (see [[measure_gt]])""")

c.add_property(name = "F_tau",
               getter = cfunction("block_gf_view<triqs::mesh::imtime> F_tau ()"),
               doc = r"""Improved estimator function $F^\sigma_{ab}(\tau)$ (see [[measure_gt]])""")

c.add_property(name = "nn_tau",
               getter = cfunction("block_gf_view<triqs::mesh::imtime> nn_tau ()"),
               doc = r"""Density-density correlation function $\langle n^\sigma_{a}(\tau)
   n^{\sigma'}_{b}(0)\rangle$ (see [[measure_nnt]])""")

c.add_property(name = "G_iw",
               getter = cfunction("block_gf_view<triqs::mesh::imfreq> G_iw ()"),
               doc = r"""Impurity Green's function $G^\sigma_{ab}(i\omega)$ (see [[measure_gw]])""")

c.add_property(name = "F_iw",
               getter = cfunction("block_gf_view<triqs::mesh::imfreq> F_iw ()"),
               doc = r"""Improved estimator function $F^\sigma_{ab}(i\omega)$ (see [[measure_gw]])""")

c.add_property(name = "Sigma_iw",
               getter = cfunction("block_gf_view<triqs::mesh::imfreq> Sigma_iw ()"),
               doc = r"""Impurity self-energy $\Sigma^\sigma_{ab}(i\omega)$ (see [[measure_gw]])""")

c.add_property(name = "nn_iw",
               getter = cfunction("block_gf_view<triqs::mesh::imfreq> nn_iw ()"),
               doc = r"""Density-density correlation function $\mathrm{FT}\left[\langle
   n^\sigma_{a}(\tau) n^{\sigma'}_{b}(0)\rangle\right]$ (see [[measure_nnw]])""")

c.add_property(name = "G_2w",
               getter = cfunction("triqs_ctseg::block_f_om_nu_tv_vt G_2w ()"),
               doc = r"""3-point correlation function $\chi^{\sigma\sigma'}_{abc}(i\omega,i\Omega)$
   (see [[measure_g2w]])""")

c.add_property(name = "F_2w",
               getter = cfunction("triqs_ctseg::block_f_om_nu_tv_vt F_2w ()"),
               doc = r"""3-point improved estimator (see [[measure_g2w]])""")

c.add_property(name = "G_3w",
               getter = cfunction("triqs_ctseg::gf_3w_container_t G_3w ()"),
               doc = r"""4-point correlation function
   $\chi^{\sigma\sigma'}_{abcd}(i\omega,i\omega',i\Omega)$ (see
   [[measure_g3w]])""")

c.add_property(name = "F_3w",
               getter = cfunction("triqs_ctseg::gf_3w_container_t F_3w ()"),
               doc = r"""4-point improved estimator (see [[measure_g3w]])""")

c.add_property(name = "G_l",
               getter = cfunction("triqs_ctseg::g_l_vt G_l ()"),
               doc = r"""Impurity Green's function in Legendre basis $G^\sigma_{ab}(n)$ (see
   [[measure_gl]])""")

c.add_property(name = "F_l",
               getter = cfunction("triqs_ctseg::g_l_vt F_l ()"),
               doc = r"""Improved estimator function in Legendre basis $G^\sigma_{ab}(n)$ (see
   [[measure_gl]])""")

c.add_property(name = "chipm_tau",
               getter = cfunction("gf_view<triqs::mesh::imtime> chipm_tau ()"),
               doc = r"""Spin spin correlation function $\langle s_+(\tau) s_-(0)\rangle$ (see
   [[measure_chipmt]])""")

c.add_property(name = "nn",
               getter = cfunction("block_matrix<double> nn ()"),
               doc = r"""density-density static correlation $\langle n^\sigma_a n^{\sigma'}_b
   \rangle$ (see [[measure_nn]])""")

c.add_property(name = "histogram",
               getter = cfunction("matrix_view<double> histogram ()"),
               doc = r"""histogram of hybridization perturbation order (see [[measure_hist]])""")

c.add_property(name = "histogram_composite",
               getter = cfunction("matrix_view<double> histogram_composite ()"),
               doc = r"""histogram of $\mathcal{J}_\perp$ perturbation order (see
   [[measure_hist_composite]])""")

c.add_property(name = "state_histogram",
               getter = cfunction("vector_view<double> state_histogram ()"),
               doc = r"""histogram of the boundary states of the trace (see [[measure_statehist]])""")

c.add_property(name = "average_sign",
               getter = cfunction("double average_sign ()"),
               doc = r"""Monte Carlo sign""")

c.add_property(name = "percent_done",
               getter = cfunction("double percent_done ()"),
               doc = r"""""")

module.add_class(c)


# Converter for solve_params_t
c = converter_(
        c_type = "triqs_ctseg::solve_params_t",
        doc = r"""""",
)
c.add_member(c_name = "h_int",
             c_type = "triqs_ctseg::Op",
             initializer = """  """,
             doc = r"""local Hamiltonian""")

c.add_member(c_name = "n_cycles",
             c_type = "int",
             initializer = """  """,
             doc = r"""Number of QMC cycles""")

c.add_member(c_name = "length_cycle",
             c_type = "int",
             initializer = """ 50 """,
             doc = r"""Length of a single QMC cycle""")

c.add_member(c_name = "n_warmup_cycles",
             c_type = "int",
             initializer = """ 5000 """,
             doc = r"""Number of cycles for thermalization""")

c.add_member(c_name = "random_seed",
             c_type = "int",
             initializer = """ 34788+928374*mpi::communicator().rank() """,
             doc = r"""Seed for random number generator""")

c.add_member(c_name = "random_name",
             c_type = "std::string",
             initializer = """ "" """,
             doc = r"""Name of random number generator""")

c.add_member(c_name = "max_time",
             c_type = "int",
             initializer = """ -1 """,
             doc = r"""Maximum runtime in seconds, use -1 to set infinite""")

c.add_member(c_name = "verbosity",
             c_type = "int",
             initializer = """ mpi::communicator().rank()==0?3:0 """,
             doc = r"""Verbosity level""")

c.add_member(c_name = "move_insert_segment",
             c_type = "bool",
             initializer = """ true """,
             doc = r"""Whether to perform the move insert segment (see [[move_insert_segment]])""")

c.add_member(c_name = "move_remove_segment",
             c_type = "bool",
             initializer = """ true """,
             doc = r"""Whether to perform the move remove segment (see [[move_remove_segment]])""")

c.add_member(c_name = "move_move",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Whether to perform the move move segment (see [[move_move]])""")

c.add_member(c_name = "move_swap_empty_lines",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Whether to perform the move swap empty lines (see
   [[move_swap_empty_lines]])""")

c.add_member(c_name = "move_group_into_spin_segment",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Whether to perform the move group into spin segment (see
   [[move_group_into_spin_segment]])""")

c.add_member(c_name = "move_split_spin_segment",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Whether to perform the move split spin segment (see
   [[move_split_spin_segment]])""")

c.add_member(c_name = "move_group_into_spin_segment2",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Whether to perform the move group into spin segment (see
   [[move_group_into_spin_segment2]])""")

c.add_member(c_name = "move_split_spin_segment2",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Whether to perform the move split spin segment (see
   [[move_split_spin_segment2]])""")

c.add_member(c_name = "keep_Jperp_negative",
             c_type = "bool",
             initializer = """ true """,
             doc = r"""Whether to keep Jperp negative""")

c.add_member(c_name = "measure_gt",
             c_type = "bool",
             initializer = """ true """,
             doc = r"""Whether to measure G(tau) (see [[measure_gt]])""")

c.add_member(c_name = "measure_sign",
             c_type = "bool",
             initializer = """ true """,
             doc = r"""Whether to measure MC sign (see [[measure_sign]])""")

c.add_member(c_name = "measure_ft",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Whether to measure improved estimator F(tau) (see [[measure_gt]]) (only
   isotropic spin-spin interactions if any)""")

c.add_member(c_name = "measure_gl",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Whether to measure G(l) (Legendre) (see [[measure_gl]])""")

c.add_member(c_name = "measure_fl",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Whether to measure improved estimator F(l) (Legendre) (see [[measure_gl]])
   (only isotropic spin-spin interactions if any)""")

c.add_member(c_name = "measure_gw",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Whether to measure G(iomega) (see [[measure_gw]])""")

c.add_member(c_name = "use_nfft_for_gw",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Whether to use NFFT in the measurement of gw""")

c.add_member(c_name = "measure_fw",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Whether to measure improved estimator F(iomega) (see [[measure_gw]])(only
   isotropic spin-spin interactions if any)""")

c.add_member(c_name = "measure_g2w",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Whether to measure two-frequency correlation function (see
   [[measure_g2w]])""")

c.add_member(c_name = "use_nfft_for_Mw",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Whether to use NFFT for the precomputation of Mw""")

c.add_member(c_name = "measure_f2w",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Whether to measure two-frequency improved estimator (see [[measure_g2w]])""")

c.add_member(c_name = "measure_g3w",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Whether to measure three-frequency correlation function (see
   [[measure_g3w]])""")

c.add_member(c_name = "measure_f3w",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Whether to measure three-frequency improved estimator (see
   [[measure_g3w]])""")

c.add_member(c_name = "measure_chipmt",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Whether to measure chi_{+-}(tau) (see [[measure_chipmt]])""")

c.add_member(c_name = "measure_nn",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Whether to measure <nn> (see [[measure_nn]])""")

c.add_member(c_name = "measure_nnt",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Whether to measure langle n(tau)n(0)rangle (see [[measure_nnt]])""")

c.add_member(c_name = "measure_nnw",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Whether to measure chi(iomega) (see [[measure_nnw]])""")

c.add_member(c_name = "evaluate_vertex",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Whether to evaluate vertex functions (see [[evaluate_3w_vertex]])""")

c.add_member(c_name = "measure_hist",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Whether to measure the perturbation order histogram (see [[measure_hist]])""")

c.add_member(c_name = "measure_hist_composite",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Whether to measure the perturbation order histogram for bosonic lines (see
   [[measure_hist_composite]])""")

c.add_member(c_name = "measure_statehist",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Whether to measure histogram of impurity states (see
   [[measure_statehist]])""")

c.add_member(c_name = "n_w_f_vertex",
             c_type = "int",
             initializer = """ 10 """,
             doc = r"""Number of fermionic Matsubara frequencies for vertex functions""")

c.add_member(c_name = "n_w_b_vertex",
             c_type = "int",
             initializer = """ 10 """,
             doc = r"""Number of bosonic Matsubara frequencies for vertex functions""")

c.add_member(c_name = "fname_gammaw",
             c_type = "std::string",
             initializer = """ "gammaw.h5" """,
             doc = r"""File name for 4-leg vertex gammaw""")

c.add_member(c_name = "nfft_threshold",
             c_type = "int",
             initializer = """ 0 """,
             doc = r"""not to have scripts fail if code is compiled without NFFT support
   Warning: parameter nfft_threshold will not be used because the code has
   been compiled without NFFT support""")

c.add_member(c_name = "hartree_shift",
             c_type = "std::vector<double>",
             initializer = """ std::vector<double>{} """,
             doc = r"""Hartree shift of the chem pot""")

module.add_converter(c)

# Converter for constr_params_t
c = converter_(
        c_type = "triqs_ctseg::constr_params_t",
        doc = r"""""",
)
c.add_member(c_name = "beta",
             c_type = "double",
             initializer = """  """,
             doc = r"""Inverse temperature""")

c.add_member(c_name = "gf_struct",
             c_type = "triqs_ctseg::gf_struct_t",
             initializer = """  """,
             doc = r"""Structure of the GF (names, sizes of blocks)""")

c.add_member(c_name = "n_tau",
             c_type = "int",
             initializer = """ 10001 """,
             doc = r"""Number of time slices for $Delta(\tau)$/$G(\tau)$/$F(\tau)$""")

c.add_member(c_name = "n_tau_k",
             c_type = "int",
             initializer = """ 10001 """,
             doc = r"""Number of time slices for $K(\tau)$""")

c.add_member(c_name = "n_tau_jperp",
             c_type = "int",
             initializer = """ 10001 """,
             doc = r"""Number of time slices for $J_\perp(\tau)$""")

c.add_member(c_name = "n_tau_nn",
             c_type = "int",
             initializer = """ 101 """,
             doc = r"""Number of Legendre coefficients for G(l)""")

c.add_member(c_name = "n_w_b_nn",
             c_type = "int",
             initializer = """ 32 """,
             doc = r"""Number of bosonic Matsubara frequencies for $nn(i\omega)$,
   $\mathcal{D}_0(i\omega)$, $J_\perp(i\omega)$""")

c.add_member(c_name = "n_iw",
             c_type = "int",
             initializer = """ 1025 """,
             doc = r"""Number of fermionic Matsubara frequencies for $G_0(i\omega)$, $G$, $F$,
   $\Sigma$""")

c.add_member(c_name = "n_legendre_g",
             c_type = "int",
             initializer = """ 256 """,
             doc = r"""""")

module.add_converter(c)


module.generate_code()
