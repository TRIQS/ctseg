
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

#define C2PY_VERSION_MAJOR 1
#define C2PY_VERSION_MINOR 0

#include <c2py/c2py.hpp>
#include <c2py/version_check.hpp>
#include <c2py/serialization/h5.hpp>

using c2py::operator""_a;

// ==================== enums =====================

// ==================== module classes =====================

// --------- class _c2py_cls_84421402 -----------
using _c2py_cls_84421402                                            = triqs_ctseg::constr_params_t;
template <> constexpr bool c2py::is_wrapped<_c2py_cls_84421402>     = true;
template <> inline constexpr auto c2py::tp_name<_c2py_cls_84421402> = "triqs_ctseg.solver_core.ConstrParamsT";

static int synth_constructor_a9a3d741(PyObject *self, PyObject *args, PyObject *kwargs) {
  if (args and PyTuple_Check(args) and (PyTuple_Size(args) > 0)) {
    PyErr_SetString(
       PyExc_RuntimeError,
       ("Error in constructing triqs_ctseg::constr_params_t.\nNo positional arguments allowed. Use keywords arguments"));
    return -1;
  }
  c2py::pydict_extractor de{kwargs};
  try {
    ((c2py::wrap<_c2py_cls_84421402> *)self)->_c = new _c2py_cls_84421402{};
  } catch (std::exception const &e) {
    PyErr_SetString(
       PyExc_RuntimeError,
       ("Error in constructing triqs_ctseg::constr_params_t from a Python dict.\n   "s + e.what()).c_str());
    return -1;
  }
  auto &self_c = *(((c2py::wrap<_c2py_cls_84421402> *)self)->_c);
  de("beta", self_c.beta, false);
  de("gf_struct", self_c.gf_struct, false);
  de("n_tau", self_c.n_tau, true);
  de("n_tau_bosonic", self_c.n_tau_bosonic, true);
  return de.check();
}

template <> constexpr initproc c2py::tp_init<_c2py_cls_84421402> = synth_constructor_a9a3d741;

template <>
const std::string c2py::tp_ctor_doc<_c2py_cls_84421402> =
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
// clang-format off
template <> PyMethodDef c2py::tp_methods<_c2py_cls_84421402>[] = {
   {nullptr, nullptr, 0, nullptr} // Sentinel
};
// clang-format on

constexpr auto _c2py_doc_member_c1c60074 = R"DOC(Inverse temperature :math:`\beta`.)DOC";
constexpr auto _c2py_doc_member_0848bcf9 = R"DOC(Structure of the Green's function (names and sizes of blocks).)DOC";
constexpr auto _c2py_doc_member_e8829399 = R"DOC(Number of time slices for fermionic functions.)DOC";
constexpr auto _c2py_doc_member_95669fd7 = R"DOC(Number of time slices for bosonic functions.)DOC";
static PyObject *prop_get_dict_a9a3d741(PyObject *self, void *) {
  auto &self_c = *(((c2py::wrap<_c2py_cls_84421402> *)self)->_c);
  c2py::pydict dic;
  dic["beta"]          = self_c.beta;
  dic["gf_struct"]     = self_c.gf_struct;
  dic["n_tau"]         = self_c.n_tau;
  dic["n_tau_bosonic"] = self_c.n_tau_bosonic;
  return dic.new_ref();
}

// ----- Member and property table ----

template <>
constinit PyGetSetDef c2py::tp_getset<_c2py_cls_84421402>[] = {
   c2py::getsetdef_from_member<&_c2py_cls_84421402::beta, _c2py_cls_84421402>("beta", _c2py_doc_member_c1c60074),
   c2py::getsetdef_from_member<&_c2py_cls_84421402::gf_struct, _c2py_cls_84421402>("gf_struct",
                                                                                   _c2py_doc_member_0848bcf9),
   c2py::getsetdef_from_member<&_c2py_cls_84421402::n_tau, _c2py_cls_84421402>("n_tau", _c2py_doc_member_e8829399),
   c2py::getsetdef_from_member<&_c2py_cls_84421402::n_tau_bosonic, _c2py_cls_84421402>("n_tau_bosonic",
                                                                                       _c2py_doc_member_95669fd7),
   {"__dict__", (getter)prop_get_dict_a9a3d741, nullptr, "", nullptr},
   {nullptr, nullptr, nullptr, nullptr, nullptr}};

template <>
const std::string c2py::tp_doc<_c2py_cls_84421402> = R"DOC(Parameters used for constructing the solver class.)DOC"
   + std::string{"\n\n----------\n\n"} + c2py::tp_ctor_doc<_c2py_cls_84421402>;
// --------- class _c2py_cls_854eef0c -----------
using _c2py_cls_854eef0c                                            = triqs_ctseg::solve_params_t;
template <> constexpr bool c2py::is_wrapped<_c2py_cls_854eef0c>     = true;
template <> inline constexpr auto c2py::tp_name<_c2py_cls_854eef0c> = "triqs_ctseg.solver_core.SolveParamsT";

static int synth_constructor_9708850a(PyObject *self, PyObject *args, PyObject *kwargs) {
  if (args and PyTuple_Check(args) and (PyTuple_Size(args) > 0)) {
    PyErr_SetString(
       PyExc_RuntimeError,
       ("Error in constructing triqs_ctseg::solve_params_t.\nNo positional arguments allowed. Use keywords arguments"));
    return -1;
  }
  c2py::pydict_extractor de{kwargs};
  try {
    ((c2py::wrap<_c2py_cls_854eef0c> *)self)->_c = new _c2py_cls_854eef0c{};
  } catch (std::exception const &e) {
    PyErr_SetString(PyExc_RuntimeError,
                    ("Error in constructing triqs_ctseg::solve_params_t from a Python dict.\n   "s + e.what()).c_str());
    return -1;
  }
  auto &self_c = *(((c2py::wrap<_c2py_cls_854eef0c> *)self)->_c);
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
  de("length_cycle_g3", self_c.length_cycle_g3, true);
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
  de("measure_state_hist", self_c.measure_state_hist, true);
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

template <> constexpr initproc c2py::tp_init<_c2py_cls_854eef0c> = synth_constructor_9708850a;

template <>
const std::string c2py::tp_ctor_doc<_c2py_cls_854eef0c> =
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

length_cycle_g3 : {par_10}, default=0

n_warmup_cycles : {par_11}, default=5000

random_seed : {par_12}, default=34788 + 928374 * mpi::communicator().rank()

random_name : {par_13}, default=""

max_time : {par_14}, default=-1

verbosity : {par_15}, default== 0 ? 3 : 0

move_insert_segment : {par_16}, default=true

move_remove_segment : {par_17}, default=true

move_double_insert_segment : {par_18}, default=true

move_double_remove_segment : {par_19}, default=true

move_move_segment : {par_20}, default=true

move_split_segment : {par_21}, default=true

move_regroup_segment : {par_22}, default=true

move_insert_spin_segment : {par_23}, default=true

move_remove_spin_segment : {par_24}, default=true

move_split_spin_segment : {par_25}, default=true

move_regroup_spin_segment : {par_26}, default=true

move_swap_spin_lines : {par_27}, default=true

measure_pert_order : {par_28}, default=true

measure_G_tau : {par_29}, default=true

measure_F_tau : {par_30}, default=false

measure_densities : {par_31}, default=true

measure_average_sign : {par_32}, default=true

measure_nn_static : {par_33}, default=false

measure_nn_tau : {par_34}, default=false

measure_nn_nu_dlr : {par_35}, default=false

measure_Sperp_tau : {par_36}, default=false

measure_state_hist : {par_37}, default=false

measure_g2w : {par_38}, default=false

measure_g3w : {par_39}, default=false

imag_threshold : {par_40}, default=1.e-13

det_init_size : {par_41}, default=100

det_n_operations_before_check : {par_42}, default=100

det_precision_warning : {par_43}, default=1.e-8

det_precision_error : {par_44}, default=1.e-5

det_singular_threshold : {par_45}, default=-1

histogram_max_order : {par_46}, default=1000

visualize_config : {par_47}, default=false

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
                       c2py::python_typename<bool>()});

// ----- Method table ----
// clang-format off
template <> PyMethodDef c2py::tp_methods<_c2py_cls_854eef0c>[] = {
   {nullptr, nullptr, 0, nullptr} // Sentinel
};
// clang-format on

constexpr auto _c2py_doc_member_d8aa2356 = R"DOC(Quartic part of the local Hamiltonian.)DOC";
constexpr auto _c2py_doc_member_ecdc267f =
   R"DOC(Quandratic part of the local Hamiltonian (including chemical potential).)DOC";
constexpr auto _c2py_doc_member_ed24eebd =
   R"DOC(Number of points on which to measure :math:`G(\tau)` / :math:`F(\tau)` (defaults to ``n_tau``).)DOC";
constexpr auto _c2py_doc_member_4f09f16c =
   R"DOC(Number of points on which to measure 2-point functions (defaults to ``n_tau_bosonic``.))DOC";
constexpr auto _c2py_doc_member_49f82fff = R"DOC(DLR frequency cutoff.)DOC";
constexpr auto _c2py_doc_member_8827d563 = R"DOC(DLR precision.)DOC";
constexpr auto _c2py_doc_member_fb26f8b4 =
   R"DOC(Number of bosonic M-frequency points on which to measure vertex functions.)DOC";
constexpr auto _c2py_doc_member_807e0f18 =
   R"DOC(Number of fermionic M-frequency points on which to measure vertex functions.)DOC";
constexpr auto _c2py_doc_member_14be584c = R"DOC(Number of QMC cycles.)DOC";
constexpr auto _c2py_doc_member_70b576e5 = R"DOC(Length of a single QMC cycle.)DOC";
constexpr auto _c2py_doc_member_4cee6c28 =
   R"DOC(Number of QMC moves between two measurements of g3w (defaults to ``length_cycle``).
Must be a multiple of ``length_cycle``.)DOC";
constexpr auto _c2py_doc_member_fda38c0f = R"DOC(Number of cycles for thermalization.)DOC";
constexpr auto _c2py_doc_member_f20d74b5 = R"DOC(Seed for random number generator.)DOC";
constexpr auto _c2py_doc_member_36d8a6c3 = R"DOC(Name of random number generator.)DOC";
constexpr auto _c2py_doc_member_d3318c92 = R"DOC(Maximum runtime in seconds, use -1 to set infinite.)DOC";
constexpr auto _c2py_doc_member_60f22021 = R"DOC(Verbosity level.)DOC";
constexpr auto _c2py_doc_member_18467b51 = R"DOC(Whether to perform the move insert segment.)DOC";
constexpr auto _c2py_doc_member_fee2a3c0 = R"DOC(Whether to perform the move remove segment.)DOC";
constexpr auto _c2py_doc_member_fb0e39cb = R"DOC(Whether to perform the move double insert segment.)DOC";
constexpr auto _c2py_doc_member_14f52b3a = R"DOC(Whether to perform the move double remove segment.)DOC";
constexpr auto _c2py_doc_member_d2971765 = R"DOC(Whether to perform the move move segment.)DOC";
constexpr auto _c2py_doc_member_e799aa44 = R"DOC(Whether to perform the move split segment.)DOC";
constexpr auto _c2py_doc_member_2146e9f0 = R"DOC(Whether to perform the move group into spin segment.)DOC";
constexpr auto _c2py_doc_member_089e9486 = R"DOC(Whether to perform the move insert spin segment.)DOC";
constexpr auto _c2py_doc_member_227cd59d = R"DOC(Whether to perform the move remove spin segment.)DOC";
constexpr auto _c2py_doc_member_37bb2099 = R"DOC(Whether to perform the move insert spin segment.)DOC";
constexpr auto _c2py_doc_member_001fd00d = R"DOC(Whether to perform the move remove spin segment.)DOC";
constexpr auto _c2py_doc_member_9bd3e06a = R"DOC(Whether to perform the move swap spin lines.)DOC";
constexpr auto _c2py_doc_member_abed454f =
   R"DOC(Whether to measure the perturbation order histograms (order in Delta and Jperp).)DOC";
constexpr auto _c2py_doc_member_5c2ce6b3 = R"DOC(Whether to measure :math:`G(\tau)`.)DOC";
constexpr auto _c2py_doc_member_b49a0194 = R"DOC(Whether to measure :math:`F(\tau)`.)DOC";
constexpr auto _c2py_doc_member_a21447db = R"DOC(Whether to measure densities.)DOC";
constexpr auto _c2py_doc_member_01ec2634 = R"DOC(Whether to measure the average sign.)DOC";
constexpr auto _c2py_doc_member_ba037394 = R"DOC(Whether to measure :math:`\langle n(0) n(0) \rangle`.)DOC";
constexpr auto _c2py_doc_member_e0376650 = R"DOC(Whether to measure :math:`\langle n(\tau) n(0) \rangle`.)DOC";
constexpr auto _c2py_doc_member_3870093e = R"DOC(Whether to measure :math:`\langle n(\nu)n(0) \rangle`.)DOC";
constexpr auto _c2py_doc_member_e210e472 = R"DOC(Whether to measure :math:`\langle S_x(\tau) S_x(0) \rangle`.)DOC";
constexpr auto _c2py_doc_member_2d047c77 = R"DOC(Whether to measure state histograms.)DOC";
constexpr auto _c2py_doc_member_58565711 = R"DOC(Whether to measure three-point correlation function.)DOC";
constexpr auto _c2py_doc_member_585895a8 = R"DOC(Whether to measure four-point correlation function.)DOC";
constexpr auto _c2py_doc_member_74ea2d22 =
   R"DOC(Threshold below which the imaginary part of the local Hamiltonian h_loc0 is set to zero
(CT-SEG uses a real h_loc0); above it the solver errors. Raise to accept a larger
imaginary part.)DOC";
constexpr auto _c2py_doc_member_c5a59be0 = R"DOC(The maximum size of the determinant matrix before a resize.)DOC";
constexpr auto _c2py_doc_member_57c4615c =
   R"DOC(Max number of ops before testing the accuracy of :math:`\det(M)` and :math:`M^{-1}`.)DOC";
constexpr auto _c2py_doc_member_ea92d1b1 = R"DOC(Threshold for determinant precision warnings.)DOC";
constexpr auto _c2py_doc_member_374c1d87 = R"DOC(Threshold for determinant precision error.)DOC";
constexpr auto _c2py_doc_member_2e4d5157 =
   R"DOC(Bound for the determinant matrix being singular (if :math:`< 0`, checks for subnormal numbers).)DOC";
constexpr auto _c2py_doc_member_08ca4150 = R"DOC(Maximum order for the perturbation order histograms.)DOC";
constexpr auto _c2py_doc_member_1c530ed9 = R"DOC(Output characteristic configurations in a separate file.)DOC";
static PyObject *prop_get_dict_9708850a(PyObject *self, void *) {
  auto &self_c = *(((c2py::wrap<_c2py_cls_854eef0c> *)self)->_c);
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
  dic["length_cycle_g3"]               = self_c.length_cycle_g3;
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
  dic["measure_state_hist"]            = self_c.measure_state_hist;
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
constinit PyGetSetDef c2py::tp_getset<_c2py_cls_854eef0c>[] = {
   c2py::getsetdef_from_member<&_c2py_cls_854eef0c::h_int, _c2py_cls_854eef0c>("h_int", _c2py_doc_member_d8aa2356),
   c2py::getsetdef_from_member<&_c2py_cls_854eef0c::h_loc0, _c2py_cls_854eef0c>("h_loc0", _c2py_doc_member_ecdc267f),
   c2py::getsetdef_from_member<&_c2py_cls_854eef0c::n_tau_G, _c2py_cls_854eef0c>("n_tau_G", _c2py_doc_member_ed24eebd),
   c2py::getsetdef_from_member<&_c2py_cls_854eef0c::n_tau_chi2, _c2py_cls_854eef0c>("n_tau_chi2",
                                                                                    _c2py_doc_member_4f09f16c),
   c2py::getsetdef_from_member<&_c2py_cls_854eef0c::dlr_omega_max, _c2py_cls_854eef0c>("dlr_omega_max",
                                                                                       _c2py_doc_member_49f82fff),
   c2py::getsetdef_from_member<&_c2py_cls_854eef0c::dlr_epsilon, _c2py_cls_854eef0c>("dlr_epsilon",
                                                                                     _c2py_doc_member_8827d563),
   c2py::getsetdef_from_member<&_c2py_cls_854eef0c::n_w_b_vertex, _c2py_cls_854eef0c>("n_w_b_vertex",
                                                                                      _c2py_doc_member_fb26f8b4),
   c2py::getsetdef_from_member<&_c2py_cls_854eef0c::n_w_f_vertex, _c2py_cls_854eef0c>("n_w_f_vertex",
                                                                                      _c2py_doc_member_807e0f18),
   c2py::getsetdef_from_member<&_c2py_cls_854eef0c::n_cycles, _c2py_cls_854eef0c>("n_cycles",
                                                                                  _c2py_doc_member_14be584c),
   c2py::getsetdef_from_member<&_c2py_cls_854eef0c::length_cycle, _c2py_cls_854eef0c>("length_cycle",
                                                                                      _c2py_doc_member_70b576e5),
   c2py::getsetdef_from_member<&_c2py_cls_854eef0c::length_cycle_g3, _c2py_cls_854eef0c>("length_cycle_g3",
                                                                                         _c2py_doc_member_4cee6c28),
   c2py::getsetdef_from_member<&_c2py_cls_854eef0c::n_warmup_cycles, _c2py_cls_854eef0c>("n_warmup_cycles",
                                                                                         _c2py_doc_member_fda38c0f),
   c2py::getsetdef_from_member<&_c2py_cls_854eef0c::random_seed, _c2py_cls_854eef0c>("random_seed",
                                                                                     _c2py_doc_member_f20d74b5),
   c2py::getsetdef_from_member<&_c2py_cls_854eef0c::random_name, _c2py_cls_854eef0c>("random_name",
                                                                                     _c2py_doc_member_36d8a6c3),
   c2py::getsetdef_from_member<&_c2py_cls_854eef0c::max_time, _c2py_cls_854eef0c>("max_time",
                                                                                  _c2py_doc_member_d3318c92),
   c2py::getsetdef_from_member<&_c2py_cls_854eef0c::verbosity, _c2py_cls_854eef0c>("verbosity",
                                                                                   _c2py_doc_member_60f22021),
   c2py::getsetdef_from_member<&_c2py_cls_854eef0c::move_insert_segment, _c2py_cls_854eef0c>("move_insert_segment",
                                                                                             _c2py_doc_member_18467b51),
   c2py::getsetdef_from_member<&_c2py_cls_854eef0c::move_remove_segment, _c2py_cls_854eef0c>("move_remove_segment",
                                                                                             _c2py_doc_member_fee2a3c0),
   c2py::getsetdef_from_member<&_c2py_cls_854eef0c::move_double_insert_segment, _c2py_cls_854eef0c>(
      "move_double_insert_segment", _c2py_doc_member_fb0e39cb),
   c2py::getsetdef_from_member<&_c2py_cls_854eef0c::move_double_remove_segment, _c2py_cls_854eef0c>(
      "move_double_remove_segment", _c2py_doc_member_14f52b3a),
   c2py::getsetdef_from_member<&_c2py_cls_854eef0c::move_move_segment, _c2py_cls_854eef0c>("move_move_segment",
                                                                                           _c2py_doc_member_d2971765),
   c2py::getsetdef_from_member<&_c2py_cls_854eef0c::move_split_segment, _c2py_cls_854eef0c>("move_split_segment",
                                                                                            _c2py_doc_member_e799aa44),
   c2py::getsetdef_from_member<&_c2py_cls_854eef0c::move_regroup_segment, _c2py_cls_854eef0c>(
      "move_regroup_segment", _c2py_doc_member_2146e9f0),
   c2py::getsetdef_from_member<&_c2py_cls_854eef0c::move_insert_spin_segment, _c2py_cls_854eef0c>(
      "move_insert_spin_segment", _c2py_doc_member_089e9486),
   c2py::getsetdef_from_member<&_c2py_cls_854eef0c::move_remove_spin_segment, _c2py_cls_854eef0c>(
      "move_remove_spin_segment", _c2py_doc_member_227cd59d),
   c2py::getsetdef_from_member<&_c2py_cls_854eef0c::move_split_spin_segment, _c2py_cls_854eef0c>(
      "move_split_spin_segment", _c2py_doc_member_37bb2099),
   c2py::getsetdef_from_member<&_c2py_cls_854eef0c::move_regroup_spin_segment, _c2py_cls_854eef0c>(
      "move_regroup_spin_segment", _c2py_doc_member_001fd00d),
   c2py::getsetdef_from_member<&_c2py_cls_854eef0c::move_swap_spin_lines, _c2py_cls_854eef0c>(
      "move_swap_spin_lines", _c2py_doc_member_9bd3e06a),
   c2py::getsetdef_from_member<&_c2py_cls_854eef0c::measure_pert_order, _c2py_cls_854eef0c>("measure_pert_order",
                                                                                            _c2py_doc_member_abed454f),
   c2py::getsetdef_from_member<&_c2py_cls_854eef0c::measure_G_tau, _c2py_cls_854eef0c>("measure_G_tau",
                                                                                       _c2py_doc_member_5c2ce6b3),
   c2py::getsetdef_from_member<&_c2py_cls_854eef0c::measure_F_tau, _c2py_cls_854eef0c>("measure_F_tau",
                                                                                       _c2py_doc_member_b49a0194),
   c2py::getsetdef_from_member<&_c2py_cls_854eef0c::measure_densities, _c2py_cls_854eef0c>("measure_densities",
                                                                                           _c2py_doc_member_a21447db),
   c2py::getsetdef_from_member<&_c2py_cls_854eef0c::measure_average_sign, _c2py_cls_854eef0c>(
      "measure_average_sign", _c2py_doc_member_01ec2634),
   c2py::getsetdef_from_member<&_c2py_cls_854eef0c::measure_nn_static, _c2py_cls_854eef0c>("measure_nn_static",
                                                                                           _c2py_doc_member_ba037394),
   c2py::getsetdef_from_member<&_c2py_cls_854eef0c::measure_nn_tau, _c2py_cls_854eef0c>("measure_nn_tau",
                                                                                        _c2py_doc_member_e0376650),
   c2py::getsetdef_from_member<&_c2py_cls_854eef0c::measure_nn_nu_dlr, _c2py_cls_854eef0c>("measure_nn_nu_dlr",
                                                                                           _c2py_doc_member_3870093e),
   c2py::getsetdef_from_member<&_c2py_cls_854eef0c::measure_Sperp_tau, _c2py_cls_854eef0c>("measure_Sperp_tau",
                                                                                           _c2py_doc_member_e210e472),
   c2py::getsetdef_from_member<&_c2py_cls_854eef0c::measure_state_hist, _c2py_cls_854eef0c>("measure_state_hist",
                                                                                            _c2py_doc_member_2d047c77),
   c2py::getsetdef_from_member<&_c2py_cls_854eef0c::measure_g2w, _c2py_cls_854eef0c>("measure_g2w",
                                                                                     _c2py_doc_member_58565711),
   c2py::getsetdef_from_member<&_c2py_cls_854eef0c::measure_g3w, _c2py_cls_854eef0c>("measure_g3w",
                                                                                     _c2py_doc_member_585895a8),
   c2py::getsetdef_from_member<&_c2py_cls_854eef0c::imag_threshold, _c2py_cls_854eef0c>("imag_threshold",
                                                                                        _c2py_doc_member_74ea2d22),
   c2py::getsetdef_from_member<&_c2py_cls_854eef0c::det_init_size, _c2py_cls_854eef0c>("det_init_size",
                                                                                       _c2py_doc_member_c5a59be0),
   c2py::getsetdef_from_member<&_c2py_cls_854eef0c::det_n_operations_before_check, _c2py_cls_854eef0c>(
      "det_n_operations_before_check", _c2py_doc_member_57c4615c),
   c2py::getsetdef_from_member<&_c2py_cls_854eef0c::det_precision_warning, _c2py_cls_854eef0c>(
      "det_precision_warning", _c2py_doc_member_ea92d1b1),
   c2py::getsetdef_from_member<&_c2py_cls_854eef0c::det_precision_error, _c2py_cls_854eef0c>("det_precision_error",
                                                                                             _c2py_doc_member_374c1d87),
   c2py::getsetdef_from_member<&_c2py_cls_854eef0c::det_singular_threshold, _c2py_cls_854eef0c>(
      "det_singular_threshold", _c2py_doc_member_2e4d5157),
   c2py::getsetdef_from_member<&_c2py_cls_854eef0c::histogram_max_order, _c2py_cls_854eef0c>("histogram_max_order",
                                                                                             _c2py_doc_member_08ca4150),
   c2py::getsetdef_from_member<&_c2py_cls_854eef0c::visualize_config, _c2py_cls_854eef0c>("visualize_config",
                                                                                          _c2py_doc_member_1c530ed9),
   {"__dict__", (getter)prop_get_dict_9708850a, nullptr, "", nullptr},
   {nullptr, nullptr, nullptr, nullptr, nullptr}};

template <>
const std::string c2py::tp_doc<_c2py_cls_854eef0c> =
   R"DOC(Parameters passed to the ``solve()`` method of the solver class.)DOC" + std::string{"\n\n----------\n\n"}
   + c2py::tp_ctor_doc<_c2py_cls_854eef0c>;
// --------- class _c2py_cls_361be08a -----------
using _c2py_cls_361be08a                                            = triqs_ctseg::results_t;
template <> constexpr bool c2py::is_wrapped<_c2py_cls_361be08a>     = true;
template <> inline constexpr auto c2py::tp_name<_c2py_cls_361be08a> = "triqs_ctseg.solver_core.ResultsT";

static int synth_constructor_4842f0bc(PyObject *self, PyObject *args, PyObject *kwargs) {
  if (args and PyTuple_Check(args) and (PyTuple_Size(args) > 0)) {
    PyErr_SetString(
       PyExc_RuntimeError,
       ("Error in constructing triqs_ctseg::results_t.\nNo positional arguments allowed. Use keywords arguments"));
    return -1;
  }
  c2py::pydict_extractor de{kwargs};
  try {
    ((c2py::wrap<_c2py_cls_361be08a> *)self)->_c = new _c2py_cls_361be08a{};
  } catch (std::exception const &e) {
    PyErr_SetString(PyExc_RuntimeError,
                    ("Error in constructing triqs_ctseg::results_t from a Python dict.\n   "s + e.what()).c_str());
    return -1;
  }
  auto &self_c = *(((c2py::wrap<_c2py_cls_361be08a> *)self)->_c);
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
  de("g2w", self_c.g2w, false);
  de("g3w", self_c.g3w, false);
  de("average_sign", self_c.average_sign, false);
  return de.check();
}

template <> constexpr initproc c2py::tp_init<_c2py_cls_361be08a> = synth_constructor_4842f0bc;

template <>
const std::string c2py::tp_ctor_doc<_c2py_cls_361be08a> = c2py::replace_tags(
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
    c2py::python_typename<double>()});

// ----- Method table ----
// clang-format off
template <> PyMethodDef c2py::tp_methods<_c2py_cls_361be08a>[] = {
   {"__write_hdf5__", c2py::tpxx_write_h5<_c2py_cls_361be08a>, METH_VARARGS, "  "},
   {"__getstate__", c2py::getstate_h5<_c2py_cls_361be08a>, METH_NOARGS, ""},
   {"__setstate__", c2py::setstate_h5<_c2py_cls_361be08a>, METH_O, ""},
   {nullptr, nullptr, 0, nullptr} // Sentinel
};
// clang-format on

constexpr auto _c2py_doc_member_ce4be67a = R"DOC(Single-particle Green's function :math:`G(\tau)`.)DOC";
constexpr auto _c2py_doc_member_f42eee39 = R"DOC(Self-energy improved estimator :math:`F(\tau)`.)DOC";
constexpr auto _c2py_doc_member_6316ab9f =
   R"DOC(Density-density time correlation function :math:`\langle n_a(\tau) n_b(0) \rangle`.)DOC";
constexpr auto _c2py_doc_member_ee202de7 =
   R"DOC(Density-density frequency correlation function :math:`\langle n_a(i\nu) n_b(-i\nu) \rangle`.)DOC";
constexpr auto _c2py_doc_member_2d5158c7 =
   R"DOC(Perpendicular spin-spin correlation function :math:`\langle S_x(\tau) S_x(0) \rangle`.)DOC";
constexpr auto _c2py_doc_member_cecc8cd5 =
   R"DOC(Density-density static correlation function :math:`\langle n_a(0) n_b(0) \rangle`.)DOC";
constexpr auto _c2py_doc_member_ddddb62a = R"DOC(Density per color, organized by blocks.)DOC";
constexpr auto _c2py_doc_member_fe3eb8e5 = R"DOC(Delta perturbation order histogram.)DOC";
constexpr auto _c2py_doc_member_3a5db50b = R"DOC(Average Delta perturbation order.)DOC";
constexpr auto _c2py_doc_member_496ab5b0 = R"DOC(Jperp perturbation order histogram.)DOC";
constexpr auto _c2py_doc_member_8b359046 = R"DOC(Average Jperp perturbation order.)DOC";
constexpr auto _c2py_doc_member_22c7a238 = R"DOC(State histogram.)DOC";
constexpr auto _c2py_doc_member_6a60b494 = R"DOC(Three-point correlation function.)DOC";
constexpr auto _c2py_doc_member_6a5e75fd = R"DOC(Four-point correlation function.)DOC";
constexpr auto _c2py_doc_member_02f5ab2b = R"DOC(Average sign.)DOC";
static PyObject *prop_get_dict_4842f0bc(PyObject *self, void *) {
  auto &self_c = *(((c2py::wrap<_c2py_cls_361be08a> *)self)->_c);
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
  dic["g2w"]                 = self_c.g2w;
  dic["g3w"]                 = self_c.g3w;
  dic["average_sign"]        = self_c.average_sign;
  return dic.new_ref();
}

// ----- Member and property table ----

template <>
constinit PyGetSetDef c2py::tp_getset<_c2py_cls_361be08a>[] = {
   c2py::getsetdef_from_member<&_c2py_cls_361be08a::G_tau, _c2py_cls_361be08a>("G_tau", _c2py_doc_member_ce4be67a),
   c2py::getsetdef_from_member<&_c2py_cls_361be08a::F_tau, _c2py_cls_361be08a>("F_tau", _c2py_doc_member_f42eee39),
   c2py::getsetdef_from_member<&_c2py_cls_361be08a::nn_tau, _c2py_cls_361be08a>("nn_tau", _c2py_doc_member_6316ab9f),
   c2py::getsetdef_from_member<&_c2py_cls_361be08a::nn_nu_dlr, _c2py_cls_361be08a>("nn_nu_dlr",
                                                                                   _c2py_doc_member_ee202de7),
   c2py::getsetdef_from_member<&_c2py_cls_361be08a::Sperp_tau, _c2py_cls_361be08a>("Sperp_tau",
                                                                                   _c2py_doc_member_2d5158c7),
   c2py::getsetdef_from_member<&_c2py_cls_361be08a::nn_static, _c2py_cls_361be08a>("nn_static",
                                                                                   _c2py_doc_member_cecc8cd5),
   c2py::getsetdef_from_member<&_c2py_cls_361be08a::densities, _c2py_cls_361be08a>("densities",
                                                                                   _c2py_doc_member_ddddb62a),
   c2py::getsetdef_from_member<&_c2py_cls_361be08a::pert_order_Delta, _c2py_cls_361be08a>("pert_order_Delta",
                                                                                          _c2py_doc_member_fe3eb8e5),
   c2py::getsetdef_from_member<&_c2py_cls_361be08a::average_order_Delta, _c2py_cls_361be08a>("average_order_Delta",
                                                                                             _c2py_doc_member_3a5db50b),
   c2py::getsetdef_from_member<&_c2py_cls_361be08a::pert_order_Jperp, _c2py_cls_361be08a>("pert_order_Jperp",
                                                                                          _c2py_doc_member_496ab5b0),
   c2py::getsetdef_from_member<&_c2py_cls_361be08a::average_order_Jperp, _c2py_cls_361be08a>("average_order_Jperp",
                                                                                             _c2py_doc_member_8b359046),
   c2py::getsetdef_from_member<&_c2py_cls_361be08a::state_hist, _c2py_cls_361be08a>("state_hist",
                                                                                    _c2py_doc_member_22c7a238),
   c2py::getsetdef_from_member<&_c2py_cls_361be08a::g2w, _c2py_cls_361be08a>("g2w", _c2py_doc_member_6a60b494),
   c2py::getsetdef_from_member<&_c2py_cls_361be08a::g3w, _c2py_cls_361be08a>("g3w", _c2py_doc_member_6a5e75fd),
   c2py::getsetdef_from_member<&_c2py_cls_361be08a::average_sign, _c2py_cls_361be08a>("average_sign",
                                                                                      _c2py_doc_member_02f5ab2b),
   {"__dict__", (getter)prop_get_dict_4842f0bc, nullptr, "", nullptr},
   {nullptr, nullptr, nullptr, nullptr, nullptr}};

template <>
const std::string c2py::tp_doc<_c2py_cls_361be08a> =
   R"DOC(Container for all results accumulated by the CTQMC simulation.)DOC" + std::string{"\n\n----------\n\n"}
   + c2py::tp_ctor_doc<_c2py_cls_361be08a>;
// --------- class _c2py_cls_41fbbc24 -----------
using _c2py_cls_41fbbc24                                            = triqs_ctseg::solver_core;
template <> constexpr bool c2py::is_wrapped<_c2py_cls_41fbbc24>     = true;
template <> inline constexpr auto c2py::tp_name<_c2py_cls_41fbbc24> = "triqs_ctseg.solver_core.SolverCore";
static const auto _c2py_init_8b89b046 =
   c2py::dispatcher_c_kw_t{c2py::c_constructor<_c2py_cls_41fbbc24, const triqs_ctseg::constr_params_t &>("p")};
template <> constexpr initproc c2py::tp_init<_c2py_cls_41fbbc24> = c2py::pyfkw_constructor<_c2py_init_8b89b046>;
template <>
const std::string c2py::tp_ctor_doc<_c2py_cls_41fbbc24> =
   _c2py_init_8b89b046.doc(R"DOC(
Initialize the solver.

Parameters
----------
p : {par_0}
   Parameters used for constructing the solver class.
)DOC",
                           {{c2py::python_typename<const triqs_ctseg::constr_params_t &>()}});
// solve
static auto const _c2py_fun_fb019e73 = c2py::dispatcher_f_kw_t{c2py::cmethod(
   [](_c2py_cls_41fbbc24 &self, const triqs_ctseg::solve_params_t &p) -> decltype(auto) { return self.solve(p); },
   "self", "p")};

static const auto _c2py_doc_fb019e73 =
   _c2py_fun_fb019e73.doc(R"DOC(
Solve the impurity problem.

Parameters
----------
p : {par_0}
   Parameters controlling the MC simulation and measurements.
)DOC",
                          {{c2py::python_typename<const triqs_ctseg::solve_params_t &>()}});

// ----- Method table ----
// clang-format off
template <> PyMethodDef c2py::tp_methods<_c2py_cls_41fbbc24>[] = {
   PMDF("solve", fb019e73),
   {"__write_hdf5__", c2py::tpxx_write_h5<_c2py_cls_41fbbc24>, METH_VARARGS, "  "},
   {"__getstate__", c2py::getstate_h5<_c2py_cls_41fbbc24>, METH_NOARGS, ""},
   {"__setstate__", c2py::setstate_h5<_c2py_cls_41fbbc24>, METH_O, ""},
   {nullptr, nullptr, 0, nullptr} // Sentinel
};
// clang-format on

constexpr auto _c2py_doc_member_9e478182 = R"DOC(Parameters used for constructing the solver.)DOC";
constexpr auto _c2py_doc_member_60044ba2 = R"DOC(Parameters passed to the ``solve()`` method.)DOC";
constexpr auto _c2py_doc_member_9d471e9e = R"DOC(Container for all results accumulated by the CTQMC simulation.)DOC";
static constexpr auto prop_doc_81cfd2fd  = R"DOC(Dynamical density-density interaction :math:`D_0(\tau)`.)DOC";
static constexpr auto prop_doc_d261a2f3  = R"DOC(Hybridization function :math:`\Delta(\tau)`.)DOC";
static constexpr auto prop_doc_09e294be  = R"DOC(Dynamical spin-spin interaction :math:`\mathcal{J}_\perp(\tau)`.)DOC";

// ----- Member and property table ----

template <>
constinit PyGetSetDef c2py::tp_getset<_c2py_cls_41fbbc24>[] = {
   c2py::getsetdef_from_member<&_c2py_cls_41fbbc24::constr_params, _c2py_cls_41fbbc24>("constr_params",
                                                                                       _c2py_doc_member_9e478182),
   c2py::getsetdef_from_member<&_c2py_cls_41fbbc24::solve_params, _c2py_cls_41fbbc24>("solve_params",
                                                                                      _c2py_doc_member_60044ba2),
   c2py::getsetdef_from_member<&_c2py_cls_41fbbc24::results, _c2py_cls_41fbbc24>("results", _c2py_doc_member_9d471e9e),
   {"D0_tau", c2py::getter_from_method<c2py::castm<>(&triqs_ctseg::solver_core::D0_tau)>, nullptr, prop_doc_81cfd2fd,
    nullptr},
   {"Delta_tau", c2py::getter_from_method<c2py::castm<>(&triqs_ctseg::solver_core::Delta_tau)>, nullptr,
    prop_doc_d261a2f3, nullptr},
   {"Jperp_tau", c2py::getter_from_method<c2py::castm<>(&triqs_ctseg::solver_core::Jperp_tau)>, nullptr,
    prop_doc_09e294be, nullptr},
   {nullptr, nullptr, nullptr, nullptr, nullptr}};

template <>
const std::string c2py::tp_doc<_c2py_cls_41fbbc24> =
   R"DOC(Continuous-time hybridization-expansion quantum Monte Carlo solver.)DOC" + std::string{"\n\n----------\n\n"}
   + c2py::tp_ctor_doc<_c2py_cls_41fbbc24>;

// ==================== module functions ====================

//--------------------- module function table  -----------------------------

// clang-format off
static PyMethodDef module_methods[] = {
   {nullptr, nullptr, 0, nullptr}  // Sentinel
};
// clang-format on

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

  if (PyType_Ready(&c2py::wrap_pytype<_c2py_cls_84421402>) < 0) return NULL;
  if (PyType_Ready(&c2py::wrap_pytype<_c2py_cls_854eef0c>) < 0) return NULL;
  if (PyType_Ready(&c2py::wrap_pytype<_c2py_cls_361be08a>) < 0) return NULL;
  if (PyType_Ready(&c2py::wrap_pytype<_c2py_cls_41fbbc24>) < 0) return NULL;

  m = PyModule_Create(&module_def);
  if (m == NULL) return NULL;

  if (not c2py::register_internal_types()) return NULL;
#define _add_type(T, N)                                                                                                \
  if (not c2py::add_type_object_to_main<T>(N, m)) return NULL
  _add_type(_c2py_cls_84421402, "ConstrParamsT");
  _add_type(_c2py_cls_854eef0c, "SolveParamsT");
  _add_type(_c2py_cls_361be08a, "ResultsT");
  _add_type(_c2py_cls_41fbbc24, "SolverCore");
#undef _add_type

  c2py::pyref module = c2py::pyref::module("h5.formats");
  if (not module) return nullptr;
  c2py::pyref register_class = module.attr("register_class");

  register_h5_type<_c2py_cls_361be08a>(register_class);
  register_h5_type<_c2py_cls_41fbbc24>(register_class);

  return m;
}
#endif
// CLAIR_WRAP_GEN
