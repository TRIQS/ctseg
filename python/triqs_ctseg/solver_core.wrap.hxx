#include <c2py/c2py.hpp>

#ifndef C2PY_HXX_DECLARATION_solver_core_GUARDS
#define C2PY_HXX_DECLARATION_solver_core_GUARDS
template <> constexpr bool c2py::is_wrapped<triqs_ctseg::constr_params_t>     = true;
template <> inline constexpr auto c2py::tp_name<triqs_ctseg::constr_params_t> = "triqs_ctseg.solver_core.ConstrParamsT";
template <> constexpr bool c2py::is_wrapped<triqs_ctseg::solve_params_t>      = true;
template <> inline constexpr auto c2py::tp_name<triqs_ctseg::solve_params_t>  = "triqs_ctseg.solver_core.SolveParamsT";
template <> constexpr bool c2py::is_wrapped<triqs_ctseg::results_t>           = true;
template <> inline constexpr auto c2py::tp_name<triqs_ctseg::results_t>       = "triqs_ctseg.solver_core.ResultsT";
template <> constexpr bool c2py::is_wrapped<triqs_ctseg::solver_core>         = true;
template <> inline constexpr auto c2py::tp_name<triqs_ctseg::solver_core>     = "triqs_ctseg.solver_core.SolverCore";
#endif