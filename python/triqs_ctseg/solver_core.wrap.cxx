
// C.f. https://numpy.org/doc/1.21/reference/c-api/array.html#importing-the-api
#define PY_ARRAY_UNIQUE_SYMBOL _cpp2py_ARRAY_API
#ifndef CLAIR_C2PY_WRAP_GEN
#ifdef __clang__
// #pragma clang diagnostic ignored "-W#warnings"
#endif
#ifdef __GNUC__
#pragma GCC diagnostic ignored "-Wmissing-field-initializers"
#pragma GCC diagnostic ignored "-Wcast-function-type"
#pragma GCC diagnostic ignored "-Wcpp"
#endif

#define C2PY_VERSION_MAJOR 0
#define C2PY_VERSION_MINOR 1

#include <c2py/c2py.hpp>
#include <c2py/serialization/h5.hpp>

using c2py::operator""_a;

// ==================== enums =====================

// ==================== module classes =====================

// --------- class _c2py_cls_0 -----------
using _c2py_cls_0                                            = triqs_ctseg::constr_params_t;
template <> constexpr bool c2py::is_wrapped<_c2py_cls_0>     = true;
template <> inline constexpr auto c2py::tp_name<_c2py_cls_0> = "triqs_ctseg.solver_core.ConstrParamsT";

static int synth_constructor_0(PyObject *self, PyObject *args, PyObject *kwargs) {
  if (args and PyTuple_Check(args) and (PyTuple_Size(args) > 0)) {
    PyErr_SetString(
       PyExc_RuntimeError,
       ("Error in constructing triqs_ctseg::constr_params_t.\nNo positional arguments allowed. Use keywords arguments"));
    return -1;
  }
  c2py::pydict_extractor de{kwargs};
  try {
    ((c2py::wrap<_c2py_cls_0> *)self)->_c = new _c2py_cls_0{};
  } catch (std::exception const &e) {
    PyErr_SetString(
       PyExc_RuntimeError,
       ("Error in constructing triqs_ctseg::constr_params_t from a Python dict.\n   "s + e.what()).c_str());
    return -1;
  }
  auto &self_c = *(((c2py::wrap<_c2py_cls_0> *)self)->_c);
  de("beta", self_c.beta, false);
  de("gf_struct", self_c.gf_struct, false);
  de("n_tau", self_c.n_tau, true);
  de("n_tau_bosonic", self_c.n_tau_bosonic, true);
  return de.check();
}

template <> constexpr initproc c2py::tp_init<_c2py_cls_0> = synth_constructor_0;

template <>
const std::string c2py::tp_ctor_doc<_c2py_cls_0> =
   c2py::replace_tags(R"DOC(Synthesized constructor with the following keyword arguments:

Parameters
----------
beta : {par_0}

gf_struct : {par_1}

n_tau : {par_2}, default=10001

n_tau_bosonic : {par_3}, default=10001

)DOC",
                      "par",
                      {c2py::python_typename<double>(), c2py::python_typename<triqs::gfs::gf_struct_t>(),
                       c2py::python_typename<int>(), c2py::python_typename<int>()});

// ----- Method table ----
template <>
PyMethodDef c2py::tp_methods<_c2py_cls_0>[] = {

   {nullptr, nullptr, 0, nullptr} // Sentinel
};

constexpr auto _c2py_doc_member_0 = R"DOC(Inverse temperature :math:`\beta`.)DOC";
constexpr auto _c2py_doc_member_1 = R"DOC(Structure of the Green's function (names and sizes of blocks).)DOC";
constexpr auto _c2py_doc_member_2 = R"DOC(Number of time slices for fermionic functions.)DOC";
constexpr auto _c2py_doc_member_3 = R"DOC(Number of time slices for bosonic functions.)DOC";
static PyObject *prop_get_dict_0(PyObject *self, void *) {
  auto &self_c = *(((c2py::wrap<_c2py_cls_0> *)self)->_c);
  c2py::pydict dic;
  dic["beta"]          = self_c.beta;
  dic["gf_struct"]     = self_c.gf_struct;
  dic["n_tau"]         = self_c.n_tau;
  dic["n_tau_bosonic"] = self_c.n_tau_bosonic;
  return dic.new_ref();
}

// ----- Member and property table ----

template <>
constinit PyGetSetDef c2py::tp_getset<_c2py_cls_0>[] = {
   c2py::getsetdef_from_member<&_c2py_cls_0::beta, _c2py_cls_0>("beta", _c2py_doc_member_0),
   c2py::getsetdef_from_member<&_c2py_cls_0::gf_struct, _c2py_cls_0>("gf_struct", _c2py_doc_member_1),
   c2py::getsetdef_from_member<&_c2py_cls_0::n_tau, _c2py_cls_0>("n_tau", _c2py_doc_member_2),
   c2py::getsetdef_from_member<&_c2py_cls_0::n_tau_bosonic, _c2py_cls_0>("n_tau_bosonic", _c2py_doc_member_3),
   {"__dict__", (getter)prop_get_dict_0, nullptr, "", nullptr},
   {nullptr, nullptr, nullptr, nullptr, nullptr}};

template <>
const std::string c2py::tp_doc<_c2py_cls_0> = R"DOC(Parameters used for constructing the solver class.)DOC"
   + std::string{"\n\n----------\n\n"} + c2py::tp_ctor_doc<_c2py_cls_0>;
// --------- class _c2py_cls_1 -----------
using _c2py_cls_1                                            = triqs_ctseg::solve_params_t;
template <> constexpr bool c2py::is_wrapped<_c2py_cls_1>     = true;
template <> inline constexpr auto c2py::tp_name<_c2py_cls_1> = "triqs_ctseg.solver_core.SolveParamsT";

static int synth_constructor_1(PyObject *self, PyObject *args, PyObject *kwargs) {
  if (args and PyTuple_Check(args) and (PyTuple_Size(args) > 0)) {
    PyErr_SetString(
       PyExc_RuntimeError,
       ("Error in constructing triqs_ctseg::solve_params_t.\nNo positional arguments allowed. Use keywords arguments"));
    return -1;
  }
  c2py::pydict_extractor de{kwargs};
  try {
    ((c2py::wrap<_c2py_cls_1> *)self)->_c = new _c2py_cls_1{};
  } catch (std::exception const &e) {
    PyErr_SetString(PyExc_RuntimeError,
                    ("Error in constructing triqs_ctseg::solve_params_t from a Python dict.\n   "s + e.what()).c_str());
    return -1;
  }
  auto &self_c = *(((c2py::wrap<_c2py_cls_1> *)self)->_c);
  de("h_int", self_c.h_int, false);
  de("h_loc0", self_c.h_loc0, false);
  de("n_tau_G", self_c.n_tau_G, true);
  de("n_tau_chi2", self_c.n_tau_chi2, true);
  de("dlr_omega_max", self_c.dlr_omega_max, true);
  de("dlr_epsilon", self_c.dlr_epsilon, true);
  de("n_w_b_vertex", self_c.n_w_b_vertex, true);
  de("n_w_f_vertex", self_c.n_w_f_vertex, true);
  de("n_cycles", self_c.n_cycles, false);
  de("length_cycle", self_c.length_cycle, true);
  de("n_warmup_cycles", self_c.n_warmup_cycles, true);
  de("random_seed", self_c.random_seed, true);
  de("random_name", self_c.random_name, true);
  de("max_time", self_c.max_time, true);
  de("verbosity", self_c.verbosity, true);
  de("move_insert_segment", self_c.move_insert_segment, true);
  de("move_remove_segment", self_c.move_remove_segment, true);
  de("move_double_insert_segment", self_c.move_double_insert_segment, true);
  de("move_double_remove_segment", self_c.move_double_remove_segment, true);
  de("move_move_segment", self_c.move_move_segment, true);
  de("move_split_segment", self_c.move_split_segment, true);
  de("move_regroup_segment", self_c.move_regroup_segment, true);
  de("move_insert_spin_segment", self_c.move_insert_spin_segment, true);
  de("move_remove_spin_segment", self_c.move_remove_spin_segment, true);
  de("move_split_spin_segment", self_c.move_split_spin_segment, true);
  de("move_regroup_spin_segment", self_c.move_regroup_spin_segment, true);
  de("move_swap_spin_lines", self_c.move_swap_spin_lines, true);
  de("measure_pert_order", self_c.measure_pert_order, true);
  de("measure_G_tau", self_c.measure_G_tau, true);
  de("measure_F_tau", self_c.measure_F_tau, true);
  de("measure_densities", self_c.measure_densities, true);
  de("measure_average_sign", self_c.measure_average_sign, true);
  de("measure_nn_static", self_c.measure_nn_static, true);
  de("measure_nn_tau", self_c.measure_nn_tau, true);
  de("measure_nn_nu_dlr", self_c.measure_nn_nu_dlr, true);
  de("measure_Sperp_tau", self_c.measure_Sperp_tau, true);
  de("measure_density_matrix", self_c.measure_density_matrix, true);
  de("measure_state_hist", self_c.measure_state_hist, true);
  de("measure_dyn_corr", self_c.measure_dyn_corr, true);
  de("measure_g2w", self_c.measure_g2w, true);
  de("measure_g3w", self_c.measure_g3w, true);
  de("imag_threshold", self_c.imag_threshold, true);
  de("det_init_size", self_c.det_init_size, true);
  de("det_n_operations_before_check", self_c.det_n_operations_before_check, true);
  de("det_precision_warning", self_c.det_precision_warning, true);
  de("det_precision_error", self_c.det_precision_error, true);
  de("det_singular_threshold", self_c.det_singular_threshold, true);
  de("histogram_max_order", self_c.histogram_max_order, true);
  de("visualize_config", self_c.visualize_config, true);
  return de.check();
}

template <> constexpr initproc c2py::tp_init<_c2py_cls_1> = synth_constructor_1;

template <>
const std::string c2py::tp_ctor_doc<_c2py_cls_1> =
   c2py::replace_tags(R"DOC(Synthesized constructor with the following keyword arguments:

Parameters
----------
h_int : {par_0}

h_loc0 : {par_1}

n_cycles : {par_2}

n_tau_G : {par_3}, default=0

n_tau_chi2 : {par_4}, default=0

dlr_omega_max : {par_5}, default=100

dlr_epsilon : {par_6}, default=1e-8

n_w_b_vertex : {par_7}, default=10

n_w_f_vertex : {par_8}, default=10

length_cycle : {par_9}, default=50

n_warmup_cycles : {par_10}, default=5000

random_seed : {par_11}, default=34788 + 928374 * mpi::communicator().rank()

random_name : {par_12}, default=""

max_time : {par_13}, default=-1

verbosity : {par_14}, default== 0 ? 3 : 0

move_insert_segment : {par_15}, default=true

move_remove_segment : {par_16}, default=true

move_double_insert_segment : {par_17}, default=true

move_double_remove_segment : {par_18}, default=true

move_move_segment : {par_19}, default=true

move_split_segment : {par_20}, default=true

move_regroup_segment : {par_21}, default=true

move_insert_spin_segment : {par_22}, default=true

move_remove_spin_segment : {par_23}, default=true

move_split_spin_segment : {par_24}, default=true

move_regroup_spin_segment : {par_25}, default=true

move_swap_spin_lines : {par_26}, default=true

measure_pert_order : {par_27}, default=true

measure_G_tau : {par_28}, default=true

measure_F_tau : {par_29}, default=false

measure_densities : {par_30}, default=true

measure_average_sign : {par_31}, default=true

measure_nn_static : {par_32}, default=false

measure_nn_tau : {par_33}, default=false

measure_nn_nu_dlr : {par_34}, default=false

measure_Sperp_tau : {par_35}, default=false

measure_state_hist : {par_36}, default=false

measure_g2w : {par_37}, default=false

measure_g3w : {par_38}, default=false

imag_threshold : {par_39}, default=1.e-13

det_init_size : {par_40}, default=100

det_n_operations_before_check : {par_41}, default=100

det_precision_warning : {par_42}, default=1.e-8

det_precision_error : {par_43}, default=1.e-5

det_singular_threshold : {par_44}, default=-1

histogram_max_order : {par_45}, default=1000

visualize_config : {par_46}, default=false

measure_density_matrix : {par_47}, default=true

measure_dyn_corr : {par_48}, default=false

)DOC",
                      "par",
                      {c2py::python_typename<triqs::operators::many_body_operator>(),
                       c2py::python_typename<triqs::operators::many_body_operator>(),
                       c2py::python_typename<int>(),
                       c2py::python_typename<int>(),
                       c2py::python_typename<int>(),
                       c2py::python_typename<double>(),
                       c2py::python_typename<double>(),
                       c2py::python_typename<int>(),
                       c2py::python_typename<int>(),
                       c2py::python_typename<int>(),
                       c2py::python_typename<int>(),
                       c2py::python_typename<int>(),
                       c2py::python_typename<std::string>(),
                       c2py::python_typename<int>(),
                       c2py::python_typename<int>(),
                       c2py::python_typename<bool>(),
                       c2py::python_typename<bool>(),
                       c2py::python_typename<bool>(),
                       c2py::python_typename<bool>(),
                       c2py::python_typename<bool>(),
                       c2py::python_typename<bool>(),
                       c2py::python_typename<bool>(),
                       c2py::python_typename<bool>(),
                       c2py::python_typename<bool>(),
                       c2py::python_typename<bool>(),
                       c2py::python_typename<bool>(),
                       c2py::python_typename<bool>(),
                       c2py::python_typename<bool>(),
                       c2py::python_typename<bool>(),
                       c2py::python_typename<bool>(),
                       c2py::python_typename<bool>(),
                       c2py::python_typename<bool>(),
                       c2py::python_typename<bool>(),
                       c2py::python_typename<bool>(),
                       c2py::python_typename<bool>(),
                       c2py::python_typename<bool>(),
                       c2py::python_typename<bool>(),
                       c2py::python_typename<bool>(),
                       c2py::python_typename<bool>(),
                       c2py::python_typename<double>(),
                       c2py::python_typename<int>(),
                       c2py::python_typename<int>(),
                       c2py::python_typename<double>(),
                       c2py::python_typename<double>(),
                       c2py::python_typename<double>(),
                       c2py::python_typename<int>(),
                       c2py::python_typename<bool>(),
                       c2py::python_typename<bool>(),
                       c2py::python_typename<bool>()});

// ----- Method table ----
template <>
PyMethodDef c2py::tp_methods<_c2py_cls_1>[] = {

   {nullptr, nullptr, 0, nullptr} // Sentinel
};

constexpr auto _c2py_doc_member_4 = R"DOC(Quartic part of the local Hamiltonian.)DOC";
constexpr auto _c2py_doc_member_5 = R"DOC(Quandratic part of the local Hamiltonian (including chemical potential).)DOC";
constexpr auto _c2py_doc_member_6 =
   R"DOC(Number of points on which to measure :math:`G(\tau)` / :math:`F(\tau)` (defaults to ``n_tau``).)DOC";
constexpr auto _c2py_doc_member_7 =
   R"DOC(Number of points on which to measure 2-point functions (defaults to ``n_tau_bosonic``.))DOC";
constexpr auto _c2py_doc_member_8 = R"DOC(DLR frequency cutoff.)DOC";
constexpr auto _c2py_doc_member_9 = R"DOC(DLR precision.)DOC";
constexpr auto _c2py_doc_member_10 =
   R"DOC(Number of bosonic M-frequency points on which to measure vertex functions.)DOC";
constexpr auto _c2py_doc_member_11 =
   R"DOC(Number of fermionic M-frequency points on which to measure vertex functions.)DOC";
constexpr auto _c2py_doc_member_12 = R"DOC(Number of QMC cycles.)DOC";
constexpr auto _c2py_doc_member_13 = R"DOC(Length of a single QMC cycle.)DOC";
constexpr auto _c2py_doc_member_14 = R"DOC(Number of cycles for thermalization.)DOC";
constexpr auto _c2py_doc_member_15 = R"DOC(Seed for random number generator.)DOC";
constexpr auto _c2py_doc_member_16 = R"DOC(Name of random number generator.)DOC";
constexpr auto _c2py_doc_member_17 = R"DOC(Maximum runtime in seconds, use -1 to set infinite.)DOC";
constexpr auto _c2py_doc_member_18 = R"DOC(Verbosity level.)DOC";
constexpr auto _c2py_doc_member_19 = R"DOC(Whether to perform the move insert segment.)DOC";
constexpr auto _c2py_doc_member_20 = R"DOC(Whether to perform the move remove segment.)DOC";
constexpr auto _c2py_doc_member_21 = R"DOC(Whether to perform the move double insert segment.)DOC";
constexpr auto _c2py_doc_member_22 = R"DOC(Whether to perform the move double remove segment.)DOC";
constexpr auto _c2py_doc_member_23 = R"DOC(Whether to perform the move move segment.)DOC";
constexpr auto _c2py_doc_member_24 = R"DOC(Whether to perform the move split segment.)DOC";
constexpr auto _c2py_doc_member_25 = R"DOC(Whether to perform the move group into spin segment.)DOC";
constexpr auto _c2py_doc_member_26 = R"DOC(Whether to perform the move insert spin segment.)DOC";
constexpr auto _c2py_doc_member_27 = R"DOC(Whether to perform the move remove spin segment.)DOC";
constexpr auto _c2py_doc_member_28 = R"DOC(Whether to perform the move insert spin segment.)DOC";
constexpr auto _c2py_doc_member_29 = R"DOC(Whether to perform the move remove spin segment.)DOC";
constexpr auto _c2py_doc_member_30 = R"DOC(Whether to perform the move swap spin lines.)DOC";
constexpr auto _c2py_doc_member_31 =
   R"DOC(Whether to measure the perturbation order histograms (order in Delta and Jperp).)DOC";
constexpr auto _c2py_doc_member_32 = R"DOC(Whether to measure :math:`G(\tau)`.)DOC";
constexpr auto _c2py_doc_member_33 = R"DOC(Whether to measure :math:`F(\tau)`.)DOC";
constexpr auto _c2py_doc_member_34 = R"DOC(Whether to measure densities.)DOC";
constexpr auto _c2py_doc_member_35 = R"DOC(Whether to measure the average sign.)DOC";
constexpr auto _c2py_doc_member_36 = R"DOC(Whether to measure :math:`\langle n(0) n(0) \rangle`.)DOC";
constexpr auto _c2py_doc_member_37 = R"DOC(Whether to measure :math:`\langle n(\tau) n(0) \rangle`.)DOC";
constexpr auto _c2py_doc_member_38 = R"DOC(Whether to measure :math:`\langle n(\nu)n(0) \rangle`.)DOC";
constexpr auto _c2py_doc_member_39 = R"DOC(Whether to measure :math:`\langle S_x(\tau) S_x(0) \rangle`.)DOC";
constexpr auto _c2py_doc_member_40 = R"DOC(Whether to measure state histograms.)DOC";
constexpr auto _c2py_doc_member_41 = R"DOC(Whether to measure three-point correlation function.)DOC";
constexpr auto _c2py_doc_member_42 = R"DOC(Whether to measure four-point correlation function.)DOC";
constexpr auto _c2py_doc_member_43 =
   R"DOC(Threshold below which the imaginary part of the local Hamiltonian h_loc0 is set to zero
(CT-SEG uses a real h_loc0); above it the solver errors. Raise to accept a larger
imaginary part.)DOC";
constexpr auto _c2py_doc_member_44 = R"DOC(The maximum size of the determinant matrix before a resize.)DOC";
constexpr auto _c2py_doc_member_45 =
   R"DOC(Max number of ops before testing the accuracy of :math:`\det(M)` and :math:`M^{-1}`.)DOC";
constexpr auto _c2py_doc_member_46 = R"DOC(Threshold for determinant precision warnings.)DOC";
constexpr auto _c2py_doc_member_47 = R"DOC(Threshold for determinant precision error.)DOC";
constexpr auto _c2py_doc_member_48 =
   R"DOC(Bound for the determinant matrix being singular (if :math:`< 0`, checks for subnormal numbers).)DOC";
constexpr auto _c2py_doc_member_49 = R"DOC(Maximum order for the perturbation order histograms.)DOC";
constexpr auto _c2py_doc_member_50 = R"DOC(Output characteristic configurations in a separate file.)DOC";
constexpr auto _c2py_doc_member_70 = R"DOC(Whether to measure the occupation-basis TTI diagonal density matrix.)DOC";
constexpr auto _c2py_doc_member_71 = R"DOC(Whether to measure retarded static correlations for tail moments.)DOC";
static PyObject *prop_get_dict_1(PyObject *self, void *) {
  auto &self_c = *(((c2py::wrap<_c2py_cls_1> *)self)->_c);
  c2py::pydict dic;
  dic["h_int"]                         = self_c.h_int;
  dic["h_loc0"]                        = self_c.h_loc0;
  dic["n_tau_G"]                       = self_c.n_tau_G;
  dic["n_tau_chi2"]                    = self_c.n_tau_chi2;
  dic["dlr_omega_max"]                 = self_c.dlr_omega_max;
  dic["dlr_epsilon"]                   = self_c.dlr_epsilon;
  dic["n_w_b_vertex"]                  = self_c.n_w_b_vertex;
  dic["n_w_f_vertex"]                  = self_c.n_w_f_vertex;
  dic["n_cycles"]                      = self_c.n_cycles;
  dic["length_cycle"]                  = self_c.length_cycle;
  dic["n_warmup_cycles"]               = self_c.n_warmup_cycles;
  dic["random_seed"]                   = self_c.random_seed;
  dic["random_name"]                   = self_c.random_name;
  dic["max_time"]                      = self_c.max_time;
  dic["verbosity"]                     = self_c.verbosity;
  dic["move_insert_segment"]           = self_c.move_insert_segment;
  dic["move_remove_segment"]           = self_c.move_remove_segment;
  dic["move_double_insert_segment"]    = self_c.move_double_insert_segment;
  dic["move_double_remove_segment"]    = self_c.move_double_remove_segment;
  dic["move_move_segment"]             = self_c.move_move_segment;
  dic["move_split_segment"]            = self_c.move_split_segment;
  dic["move_regroup_segment"]          = self_c.move_regroup_segment;
  dic["move_insert_spin_segment"]      = self_c.move_insert_spin_segment;
  dic["move_remove_spin_segment"]      = self_c.move_remove_spin_segment;
  dic["move_split_spin_segment"]       = self_c.move_split_spin_segment;
  dic["move_regroup_spin_segment"]     = self_c.move_regroup_spin_segment;
  dic["move_swap_spin_lines"]          = self_c.move_swap_spin_lines;
  dic["measure_pert_order"]            = self_c.measure_pert_order;
  dic["measure_G_tau"]                 = self_c.measure_G_tau;
  dic["measure_F_tau"]                 = self_c.measure_F_tau;
  dic["measure_densities"]             = self_c.measure_densities;
  dic["measure_average_sign"]          = self_c.measure_average_sign;
  dic["measure_nn_static"]             = self_c.measure_nn_static;
  dic["measure_nn_tau"]                = self_c.measure_nn_tau;
  dic["measure_nn_nu_dlr"]             = self_c.measure_nn_nu_dlr;
  dic["measure_Sperp_tau"]             = self_c.measure_Sperp_tau;
  dic["measure_density_matrix"]        = self_c.measure_density_matrix;
  dic["measure_state_hist"]            = self_c.measure_state_hist;
  dic["measure_dyn_corr"]              = self_c.measure_dyn_corr;
  dic["measure_g2w"]                   = self_c.measure_g2w;
  dic["measure_g3w"]                   = self_c.measure_g3w;
  dic["imag_threshold"]                = self_c.imag_threshold;
  dic["det_init_size"]                 = self_c.det_init_size;
  dic["det_n_operations_before_check"] = self_c.det_n_operations_before_check;
  dic["det_precision_warning"]         = self_c.det_precision_warning;
  dic["det_precision_error"]           = self_c.det_precision_error;
  dic["det_singular_threshold"]        = self_c.det_singular_threshold;
  dic["histogram_max_order"]           = self_c.histogram_max_order;
  dic["visualize_config"]              = self_c.visualize_config;
  return dic.new_ref();
}

// ----- Member and property table ----

template <>
constinit PyGetSetDef c2py::tp_getset<_c2py_cls_1>[] = {
   c2py::getsetdef_from_member<&_c2py_cls_1::h_int, _c2py_cls_1>("h_int", _c2py_doc_member_4),
   c2py::getsetdef_from_member<&_c2py_cls_1::h_loc0, _c2py_cls_1>("h_loc0", _c2py_doc_member_5),
   c2py::getsetdef_from_member<&_c2py_cls_1::n_tau_G, _c2py_cls_1>("n_tau_G", _c2py_doc_member_6),
   c2py::getsetdef_from_member<&_c2py_cls_1::n_tau_chi2, _c2py_cls_1>("n_tau_chi2", _c2py_doc_member_7),
   c2py::getsetdef_from_member<&_c2py_cls_1::dlr_omega_max, _c2py_cls_1>("dlr_omega_max", _c2py_doc_member_8),
   c2py::getsetdef_from_member<&_c2py_cls_1::dlr_epsilon, _c2py_cls_1>("dlr_epsilon", _c2py_doc_member_9),
   c2py::getsetdef_from_member<&_c2py_cls_1::n_w_b_vertex, _c2py_cls_1>("n_w_b_vertex", _c2py_doc_member_10),
   c2py::getsetdef_from_member<&_c2py_cls_1::n_w_f_vertex, _c2py_cls_1>("n_w_f_vertex", _c2py_doc_member_11),
   c2py::getsetdef_from_member<&_c2py_cls_1::n_cycles, _c2py_cls_1>("n_cycles", _c2py_doc_member_12),
   c2py::getsetdef_from_member<&_c2py_cls_1::length_cycle, _c2py_cls_1>("length_cycle", _c2py_doc_member_13),
   c2py::getsetdef_from_member<&_c2py_cls_1::n_warmup_cycles, _c2py_cls_1>("n_warmup_cycles", _c2py_doc_member_14),
   c2py::getsetdef_from_member<&_c2py_cls_1::random_seed, _c2py_cls_1>("random_seed", _c2py_doc_member_15),
   c2py::getsetdef_from_member<&_c2py_cls_1::random_name, _c2py_cls_1>("random_name", _c2py_doc_member_16),
   c2py::getsetdef_from_member<&_c2py_cls_1::max_time, _c2py_cls_1>("max_time", _c2py_doc_member_17),
   c2py::getsetdef_from_member<&_c2py_cls_1::verbosity, _c2py_cls_1>("verbosity", _c2py_doc_member_18),
   c2py::getsetdef_from_member<&_c2py_cls_1::move_insert_segment, _c2py_cls_1>("move_insert_segment",
                                                                               _c2py_doc_member_19),
   c2py::getsetdef_from_member<&_c2py_cls_1::move_remove_segment, _c2py_cls_1>("move_remove_segment",
                                                                               _c2py_doc_member_20),
   c2py::getsetdef_from_member<&_c2py_cls_1::move_double_insert_segment, _c2py_cls_1>("move_double_insert_segment",
                                                                                      _c2py_doc_member_21),
   c2py::getsetdef_from_member<&_c2py_cls_1::move_double_remove_segment, _c2py_cls_1>("move_double_remove_segment",
                                                                                      _c2py_doc_member_22),
   c2py::getsetdef_from_member<&_c2py_cls_1::move_move_segment, _c2py_cls_1>("move_move_segment", _c2py_doc_member_23),
   c2py::getsetdef_from_member<&_c2py_cls_1::move_split_segment, _c2py_cls_1>("move_split_segment",
                                                                              _c2py_doc_member_24),
   c2py::getsetdef_from_member<&_c2py_cls_1::move_regroup_segment, _c2py_cls_1>("move_regroup_segment",
                                                                                _c2py_doc_member_25),
   c2py::getsetdef_from_member<&_c2py_cls_1::move_insert_spin_segment, _c2py_cls_1>("move_insert_spin_segment",
                                                                                    _c2py_doc_member_26),
   c2py::getsetdef_from_member<&_c2py_cls_1::move_remove_spin_segment, _c2py_cls_1>("move_remove_spin_segment",
                                                                                    _c2py_doc_member_27),
   c2py::getsetdef_from_member<&_c2py_cls_1::move_split_spin_segment, _c2py_cls_1>("move_split_spin_segment",
                                                                                   _c2py_doc_member_28),
   c2py::getsetdef_from_member<&_c2py_cls_1::move_regroup_spin_segment, _c2py_cls_1>("move_regroup_spin_segment",
                                                                                     _c2py_doc_member_29),
   c2py::getsetdef_from_member<&_c2py_cls_1::move_swap_spin_lines, _c2py_cls_1>("move_swap_spin_lines",
                                                                                _c2py_doc_member_30),
   c2py::getsetdef_from_member<&_c2py_cls_1::measure_pert_order, _c2py_cls_1>("measure_pert_order",
                                                                              _c2py_doc_member_31),
   c2py::getsetdef_from_member<&_c2py_cls_1::measure_G_tau, _c2py_cls_1>("measure_G_tau", _c2py_doc_member_32),
   c2py::getsetdef_from_member<&_c2py_cls_1::measure_F_tau, _c2py_cls_1>("measure_F_tau", _c2py_doc_member_33),
   c2py::getsetdef_from_member<&_c2py_cls_1::measure_densities, _c2py_cls_1>("measure_densities", _c2py_doc_member_34),
   c2py::getsetdef_from_member<&_c2py_cls_1::measure_average_sign, _c2py_cls_1>("measure_average_sign",
                                                                                _c2py_doc_member_35),
   c2py::getsetdef_from_member<&_c2py_cls_1::measure_nn_static, _c2py_cls_1>("measure_nn_static", _c2py_doc_member_36),
   c2py::getsetdef_from_member<&_c2py_cls_1::measure_nn_tau, _c2py_cls_1>("measure_nn_tau", _c2py_doc_member_37),
   c2py::getsetdef_from_member<&_c2py_cls_1::measure_nn_nu_dlr, _c2py_cls_1>("measure_nn_nu_dlr", _c2py_doc_member_38),
   c2py::getsetdef_from_member<&_c2py_cls_1::measure_Sperp_tau, _c2py_cls_1>("measure_Sperp_tau", _c2py_doc_member_39),
   c2py::getsetdef_from_member<&_c2py_cls_1::measure_density_matrix, _c2py_cls_1>("measure_density_matrix",
                                                                                  _c2py_doc_member_70),
   c2py::getsetdef_from_member<&_c2py_cls_1::measure_state_hist, _c2py_cls_1>("measure_state_hist",
                                                                              _c2py_doc_member_40),
   c2py::getsetdef_from_member<&_c2py_cls_1::measure_dyn_corr, _c2py_cls_1>("measure_dyn_corr",
                                                                            _c2py_doc_member_71),
   c2py::getsetdef_from_member<&_c2py_cls_1::measure_g2w, _c2py_cls_1>("measure_g2w", _c2py_doc_member_41),
   c2py::getsetdef_from_member<&_c2py_cls_1::measure_g3w, _c2py_cls_1>("measure_g3w", _c2py_doc_member_42),
   c2py::getsetdef_from_member<&_c2py_cls_1::imag_threshold, _c2py_cls_1>("imag_threshold", _c2py_doc_member_43),
   c2py::getsetdef_from_member<&_c2py_cls_1::det_init_size, _c2py_cls_1>("det_init_size", _c2py_doc_member_44),
   c2py::getsetdef_from_member<&_c2py_cls_1::det_n_operations_before_check, _c2py_cls_1>(
      "det_n_operations_before_check", _c2py_doc_member_45),
   c2py::getsetdef_from_member<&_c2py_cls_1::det_precision_warning, _c2py_cls_1>("det_precision_warning",
                                                                                 _c2py_doc_member_46),
   c2py::getsetdef_from_member<&_c2py_cls_1::det_precision_error, _c2py_cls_1>("det_precision_error",
                                                                               _c2py_doc_member_47),
   c2py::getsetdef_from_member<&_c2py_cls_1::det_singular_threshold, _c2py_cls_1>("det_singular_threshold",
                                                                                  _c2py_doc_member_48),
   c2py::getsetdef_from_member<&_c2py_cls_1::histogram_max_order, _c2py_cls_1>("histogram_max_order",
                                                                               _c2py_doc_member_49),
   c2py::getsetdef_from_member<&_c2py_cls_1::visualize_config, _c2py_cls_1>("visualize_config", _c2py_doc_member_50),
   {"__dict__", (getter)prop_get_dict_1, nullptr, "", nullptr},
   {nullptr, nullptr, nullptr, nullptr, nullptr}};

template <>
const std::string c2py::tp_doc<_c2py_cls_1> =
   R"DOC(Parameters passed to the ``solve()`` method of the solver class.)DOC" + std::string{"\n\n----------\n\n"}
   + c2py::tp_ctor_doc<_c2py_cls_1>;
// --------- class _c2py_cls_2 -----------
using _c2py_cls_2                                            = triqs_ctseg::results_t;
template <> constexpr bool c2py::is_wrapped<_c2py_cls_2>     = true;
template <> inline constexpr auto c2py::tp_name<_c2py_cls_2> = "triqs_ctseg.solver_core.ResultsT";

static int synth_constructor_2(PyObject *self, PyObject *args, PyObject *kwargs) {
  if (args and PyTuple_Check(args) and (PyTuple_Size(args) > 0)) {
    PyErr_SetString(
       PyExc_RuntimeError,
       ("Error in constructing triqs_ctseg::results_t.\nNo positional arguments allowed. Use keywords arguments"));
    return -1;
  }
  c2py::pydict_extractor de{kwargs};
  try {
    ((c2py::wrap<_c2py_cls_2> *)self)->_c = new _c2py_cls_2{};
  } catch (std::exception const &e) {
    PyErr_SetString(PyExc_RuntimeError,
                    ("Error in constructing triqs_ctseg::results_t from a Python dict.\n   "s + e.what()).c_str());
    return -1;
  }
  auto &self_c = *(((c2py::wrap<_c2py_cls_2> *)self)->_c);
  de("G_tau", self_c.G_tau, false);
  de("F_tau", self_c.F_tau, false);
  de("nn_tau", self_c.nn_tau, false);
  de("nn_nu_dlr", self_c.nn_nu_dlr, false);
  de("Sperp_tau", self_c.Sperp_tau, false);
  de("nn_static", self_c.nn_static, false);
  de("densities", self_c.densities, false);
  de("pert_order_Delta", self_c.pert_order_Delta, false);
  de("average_order_Delta", self_c.average_order_Delta, false);
  de("pert_order_Jperp", self_c.pert_order_Jperp, false);
  de("average_order_Jperp", self_c.average_order_Jperp, false);
  de("state_hist", self_c.state_hist, false);
  de("dyn_phi_n", self_c.dyn_phi_n, false);
  de("dyn_phi_phi", self_c.dyn_phi_phi, false);
  de("g2w", self_c.g2w, false);
  de("g3w", self_c.g3w, false);
  de("average_sign", self_c.average_sign, false);
  return de.check();
}

template <> constexpr initproc c2py::tp_init<_c2py_cls_2> = synth_constructor_2;

template <>
const std::string c2py::tp_ctor_doc<_c2py_cls_2> = c2py::replace_tags(
   R"DOC(Synthesized constructor with the following keyword arguments:

Parameters
----------
G_tau : {par_0}

F_tau : {par_1}

nn_tau : {par_2}

nn_nu_dlr : {par_3}

Sperp_tau : {par_4}

nn_static : {par_5}

densities : {par_6}

pert_order_Delta : {par_7}

average_order_Delta : {par_8}

pert_order_Jperp : {par_9}

average_order_Jperp : {par_10}

state_hist : {par_11}

g2w : {par_12}

g3w : {par_13}

average_sign : {par_14}

dyn_phi_n : {par_15}

dyn_phi_phi : {par_16}

)DOC",
   "par",
   {c2py::python_typename<triqs::gfs::block_gf<triqs::mesh::imtime>>(),
    c2py::python_typename<std::optional<triqs::gfs::block_gf<triqs::mesh::imtime>>>(),
    c2py::python_typename<
       std::optional<triqs::gfs::block_gf<triqs::mesh::imtime, triqs::gfs::matrix_valued, nda::C_layout, 2>>>(),
    c2py::python_typename<
       std::optional<triqs::gfs::block_gf<triqs::mesh::dlr_imfreq, triqs::gfs::matrix_valued, nda::C_layout, 2>>>(),
    c2py::python_typename<std::optional<triqs::gfs::gf<triqs::mesh::imtime>>>(),
    c2py::python_typename<std::optional<
       std::map<std::pair<std::string, std::string>,
                nda::basic_array<double, 2, nda::C_layout, 'M',
                                 nda::heap_basic<nda::mem::mallocator<nda::mem::AddressSpace::Host>>>>>>(),
    c2py::python_typename<std::optional<
       std::map<std::string,
                nda::basic_array<double, 1, nda::C_layout, 'A',
                                 nda::heap_basic<nda::mem::mallocator<nda::mem::AddressSpace::Host>>>>>>(),
    c2py::python_typename<std::optional<std::vector<double>>>(), c2py::python_typename<std::optional<double>>(),
    c2py::python_typename<std::optional<std::vector<double>>>(), c2py::python_typename<std::optional<double>>(),
    c2py::python_typename<std::optional<nda::basic_array<
       double, 1, nda::C_layout, 'V', nda::heap_basic<nda::mem::mallocator<nda::mem::AddressSpace::Host>>>>>(),
    c2py::python_typename<std::optional<triqs::gfs::block_gf<
       triqs::mesh::prod<triqs::mesh::imfreq, triqs::mesh::imfreq>, triqs::gfs::tensor_valued<4>, nda::C_layout, 2>>>(),
    c2py::python_typename<std::optional<
       triqs::gfs::block_gf<triqs::mesh::prod<triqs::mesh::imfreq, triqs::mesh::imfreq, triqs::mesh::imfreq>,
                            triqs::gfs::tensor_valued<4>, nda::C_layout, 2>>>(),
    c2py::python_typename<double>(),
    c2py::python_typename<std::optional<nda::basic_array<
       double, 2, nda::C_layout, 'M', nda::heap_basic<nda::mem::mallocator<nda::mem::AddressSpace::Host>>>>>(),
    c2py::python_typename<std::optional<nda::basic_array<
       double, 2, nda::C_layout, 'M', nda::heap_basic<nda::mem::mallocator<nda::mem::AddressSpace::Host>>>>>()});

// ----- Method table ----
template <>
PyMethodDef c2py::tp_methods<_c2py_cls_2>[] = {
   {"__write_hdf5__", c2py::tpxx_write_h5<_c2py_cls_2>, METH_VARARGS, "  "},
   {"__getstate__", c2py::getstate_h5<_c2py_cls_2>, METH_NOARGS, ""},
   {"__setstate__", c2py::setstate_h5<_c2py_cls_2>, METH_O, ""},
   {nullptr, nullptr, 0, nullptr} // Sentinel
};

constexpr auto _c2py_doc_member_51 = R"DOC(Single-particle Green's function :math:`G(\tau)`.)DOC";
constexpr auto _c2py_doc_member_52 = R"DOC(Self-energy improved estimator :math:`F(\tau)`.)DOC";
constexpr auto _c2py_doc_member_53 =
   R"DOC(Density-density time correlation function :math:`\langle n_a(\tau) n_b(0) \rangle`.)DOC";
constexpr auto _c2py_doc_member_54 =
   R"DOC(Density-density frequency correlation function :math:`\langle n_a(i\nu) n_b(-i\nu) \rangle`.)DOC";
constexpr auto _c2py_doc_member_55 =
   R"DOC(Perpendicular spin-spin correlation function :math:`\langle S_x(\tau) S_x(0) \rangle`.)DOC";
constexpr auto _c2py_doc_member_56 =
   R"DOC(Density-density static correlation function :math:`\langle n_a(0) n_b(0) \rangle`.)DOC";
constexpr auto _c2py_doc_member_57 = R"DOC(Density per color, organized by blocks.)DOC";
constexpr auto _c2py_doc_member_58 = R"DOC(Delta perturbation order histogram.)DOC";
constexpr auto _c2py_doc_member_59 = R"DOC(Average Delta perturbation order.)DOC";
constexpr auto _c2py_doc_member_60 = R"DOC(Jperp perturbation order histogram.)DOC";
constexpr auto _c2py_doc_member_61 = R"DOC(Average Jperp perturbation order.)DOC";
constexpr auto _c2py_doc_member_62 = R"DOC(State histogram.)DOC";
constexpr auto _c2py_doc_member_63 = R"DOC(Three-point correlation function.)DOC";
constexpr auto _c2py_doc_member_64 = R"DOC(Four-point correlation function.)DOC";
constexpr auto _c2py_doc_member_65 = R"DOC(Average sign.)DOC";
constexpr auto _c2py_doc_member_72 = R"DOC(Retarded source-field/density static correlation in color space.)DOC";
constexpr auto _c2py_doc_member_73 = R"DOC(Retarded source-field/source-field static correlation in color space.)DOC";
static PyObject *prop_get_dict_2(PyObject *self, void *) {
  auto &self_c = *(((c2py::wrap<_c2py_cls_2> *)self)->_c);
  c2py::pydict dic;
  dic["G_tau"]               = self_c.G_tau;
  dic["F_tau"]               = self_c.F_tau;
  dic["nn_tau"]              = self_c.nn_tau;
  dic["nn_nu_dlr"]           = self_c.nn_nu_dlr;
  dic["Sperp_tau"]           = self_c.Sperp_tau;
  dic["nn_static"]           = self_c.nn_static;
  dic["densities"]           = self_c.densities;
  dic["pert_order_Delta"]    = self_c.pert_order_Delta;
  dic["average_order_Delta"] = self_c.average_order_Delta;
  dic["pert_order_Jperp"]    = self_c.pert_order_Jperp;
  dic["average_order_Jperp"] = self_c.average_order_Jperp;
  dic["state_hist"]          = self_c.state_hist;
  dic["dyn_phi_n"]           = self_c.dyn_phi_n;
  dic["dyn_phi_phi"]         = self_c.dyn_phi_phi;
  dic["g2w"]                 = self_c.g2w;
  dic["g3w"]                 = self_c.g3w;
  dic["average_sign"]        = self_c.average_sign;
  return dic.new_ref();
}

// ----- Member and property table ----

template <>
constinit PyGetSetDef c2py::tp_getset<_c2py_cls_2>[] = {
   c2py::getsetdef_from_member<&_c2py_cls_2::G_tau, _c2py_cls_2>("G_tau", _c2py_doc_member_51),
   c2py::getsetdef_from_member<&_c2py_cls_2::F_tau, _c2py_cls_2>("F_tau", _c2py_doc_member_52),
   c2py::getsetdef_from_member<&_c2py_cls_2::nn_tau, _c2py_cls_2>("nn_tau", _c2py_doc_member_53),
   c2py::getsetdef_from_member<&_c2py_cls_2::nn_nu_dlr, _c2py_cls_2>("nn_nu_dlr", _c2py_doc_member_54),
   c2py::getsetdef_from_member<&_c2py_cls_2::Sperp_tau, _c2py_cls_2>("Sperp_tau", _c2py_doc_member_55),
   c2py::getsetdef_from_member<&_c2py_cls_2::nn_static, _c2py_cls_2>("nn_static", _c2py_doc_member_56),
   c2py::getsetdef_from_member<&_c2py_cls_2::densities, _c2py_cls_2>("densities", _c2py_doc_member_57),
   c2py::getsetdef_from_member<&_c2py_cls_2::pert_order_Delta, _c2py_cls_2>("pert_order_Delta", _c2py_doc_member_58),
   c2py::getsetdef_from_member<&_c2py_cls_2::average_order_Delta, _c2py_cls_2>("average_order_Delta",
                                                                               _c2py_doc_member_59),
   c2py::getsetdef_from_member<&_c2py_cls_2::pert_order_Jperp, _c2py_cls_2>("pert_order_Jperp", _c2py_doc_member_60),
   c2py::getsetdef_from_member<&_c2py_cls_2::average_order_Jperp, _c2py_cls_2>("average_order_Jperp",
                                                                               _c2py_doc_member_61),
   c2py::getsetdef_from_member<&_c2py_cls_2::state_hist, _c2py_cls_2>("state_hist", _c2py_doc_member_62),
   c2py::getsetdef_from_member<&_c2py_cls_2::dyn_phi_n, _c2py_cls_2>("dyn_phi_n", _c2py_doc_member_72),
   c2py::getsetdef_from_member<&_c2py_cls_2::dyn_phi_phi, _c2py_cls_2>("dyn_phi_phi", _c2py_doc_member_73),
   c2py::getsetdef_from_member<&_c2py_cls_2::g2w, _c2py_cls_2>("g2w", _c2py_doc_member_63),
   c2py::getsetdef_from_member<&_c2py_cls_2::g3w, _c2py_cls_2>("g3w", _c2py_doc_member_64),
   c2py::getsetdef_from_member<&_c2py_cls_2::average_sign, _c2py_cls_2>("average_sign", _c2py_doc_member_65),
   {"__dict__", (getter)prop_get_dict_2, nullptr, "", nullptr},
   {nullptr, nullptr, nullptr, nullptr, nullptr}};

template <>
const std::string c2py::tp_doc<_c2py_cls_2> = R"DOC(Container for all results accumulated by the CTQMC simulation.)DOC"
   + std::string{"\n\n----------\n\n"} + c2py::tp_ctor_doc<_c2py_cls_2>;
// --------- class _c2py_cls_3 -----------
using _c2py_cls_3                                            = triqs_ctseg::solver_core;
template <> constexpr bool c2py::is_wrapped<_c2py_cls_3>     = true;
template <> inline constexpr auto c2py::tp_name<_c2py_cls_3> = "triqs_ctseg.solver_core.SolverCore";
static const auto _c2py_init_0 =
   c2py::dispatcher_c_kw_t{c2py::c_constructor<_c2py_cls_3, const triqs_ctseg::constr_params_t &>("p")};
template <> constexpr initproc c2py::tp_init<_c2py_cls_3> = c2py::pyfkw_constructor<_c2py_init_0>;
template <>
const std::string c2py::tp_ctor_doc<_c2py_cls_3> =
   _c2py_init_0.doc(R"DOC(
Initialize the solver.

Parameters
----------
p : {par_0}
   Parameters used for constructing the solver class.
)DOC",
                    {{c2py::python_typename<const triqs_ctseg::constr_params_t &>()}});
// solve
static auto const _c2py_fun_0 = c2py::dispatcher_f_kw_t{c2py::cmethod(
   [](_c2py_cls_3 &self, const triqs_ctseg::solve_params_t &p) -> decltype(auto) { return self.solve(p); }, "self",
   "p")};

static const auto _c2py_doc_0 = _c2py_fun_0.doc(R"DOC(
Solve the impurity problem.

Parameters
----------
p : {par_0}
   Parameters controlling the MC simulation and measurements.
)DOC",
                                                {{c2py::python_typename<const triqs_ctseg::solve_params_t &>()}});

// ----- Method table ----
template <>
PyMethodDef c2py::tp_methods<_c2py_cls_3>[] = {
   {"solve", (PyCFunction)c2py::pyfkw<_c2py_fun_0>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_0.c_str()},
   {"__write_hdf5__", c2py::tpxx_write_h5<_c2py_cls_3>, METH_VARARGS, "  "},
   {"__getstate__", c2py::getstate_h5<_c2py_cls_3>, METH_NOARGS, ""},
   {"__setstate__", c2py::setstate_h5<_c2py_cls_3>, METH_O, ""},
   {nullptr, nullptr, 0, nullptr} // Sentinel
};

constexpr auto _c2py_doc_member_66 = R"DOC(Parameters used for constructing the solver.)DOC";
constexpr auto _c2py_doc_member_67 = R"DOC(Parameters passed to the ``solve()`` method.)DOC";
constexpr auto _c2py_doc_member_68 = R"DOC(Container for all results accumulated by the CTQMC simulation.)DOC";
static constexpr auto prop_doc_0   = R"DOC(Dynamical density-density interaction :math:`D_0(\tau)`.)DOC";
static constexpr auto prop_doc_1   = R"DOC(Hybridization function :math:`\Delta(\tau)`.)DOC";
static constexpr auto prop_doc_2   = R"DOC(Dynamical spin-spin interaction :math:`\mathcal{J}_\perp(\tau)`.)DOC";

// ----- Member and property table ----

template <>
constinit PyGetSetDef c2py::tp_getset<_c2py_cls_3>[] = {
   c2py::getsetdef_from_member<&_c2py_cls_3::constr_params, _c2py_cls_3>("constr_params", _c2py_doc_member_66),
   c2py::getsetdef_from_member<&_c2py_cls_3::solve_params, _c2py_cls_3>("solve_params", _c2py_doc_member_67),
   c2py::getsetdef_from_member<&_c2py_cls_3::results, _c2py_cls_3>("results", _c2py_doc_member_68),
   {"D0_tau", c2py::getter_from_method<c2py::castm<>(&triqs_ctseg::solver_core::D0_tau)>, nullptr, prop_doc_0, nullptr},
   {"Delta_tau", c2py::getter_from_method<c2py::castm<>(&triqs_ctseg::solver_core::Delta_tau)>, nullptr, prop_doc_1,
    nullptr},
   {"Jperp_tau", c2py::getter_from_method<c2py::castm<>(&triqs_ctseg::solver_core::Jperp_tau)>, nullptr, prop_doc_2,
    nullptr},
   {nullptr, nullptr, nullptr, nullptr, nullptr}};

template <>
const std::string c2py::tp_doc<_c2py_cls_3> =
   R"DOC(Continuous-time hybridization-expansion quantum Monte Carlo solver.)DOC" + std::string{"\n\n----------\n\n"}
   + c2py::tp_ctor_doc<_c2py_cls_3>;

// ==================== module functions ====================

//--------------------- module function table  -----------------------------

static PyMethodDef module_methods[] = {
   {nullptr, nullptr, 0, nullptr} // Sentinel
};

//--------------------- module struct & init error definition ------------

//// module doc directly in the code or "" if not present...
/// Or mandatory ?
static struct PyModuleDef module_def = {
   PyModuleDef_HEAD_INIT,
   "solver_core",                                             /* name of module */
   R"RAWDOC(Core module containing wrapped C++ code.)RAWDOC", /* module documentation, may be NULL */
   -1, /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
   module_methods,
   NULL,
   NULL,
   NULL,
   NULL};

//--------------------- module init function -----------------------------

extern "C" __attribute__((visibility("default"))) PyObject *PyInit_solver_core() {

  if (not c2py::check_python_version("solver_core")) return NULL;

  // import numpy iff 'numpy/arrayobject.h' included
#ifdef Py_ARRAYOBJECT_H
  import_array();
#endif

  PyObject *m;

  if (PyType_Ready(&c2py::wrap_pytype<c2py::py_range>) < 0) return NULL;
  if (PyType_Ready(&c2py::wrap_pytype<_c2py_cls_0>) < 0) return NULL;
  if (PyType_Ready(&c2py::wrap_pytype<_c2py_cls_1>) < 0) return NULL;
  if (PyType_Ready(&c2py::wrap_pytype<_c2py_cls_2>) < 0) return NULL;
  if (PyType_Ready(&c2py::wrap_pytype<_c2py_cls_3>) < 0) return NULL;

  m = PyModule_Create(&module_def);
  if (m == NULL) return NULL;

  auto &conv_table = *c2py::conv_table_sptr.get();

  conv_table[std::type_index(typeid(c2py::py_range)).name()] = &c2py::wrap_pytype<c2py::py_range>;
#define _add_type(T, N) c2py::add_type_object_to_main<T>(N, m, conv_table)
  _add_type(_c2py_cls_0, "ConstrParamsT");
  _add_type(_c2py_cls_1, "SolveParamsT");
  _add_type(_c2py_cls_2, "ResultsT");
  _add_type(_c2py_cls_3, "SolverCore");
#undef _add_type

  c2py::pyref module = c2py::pyref::module("h5.formats");
  if (not module) return nullptr;
  c2py::pyref register_class = module.attr("register_class");

  register_h5_type<_c2py_cls_2>(register_class);
  register_h5_type<_c2py_cls_3>(register_class);

  return m;
}
#endif
// CLAIR_WRAP_GEN
