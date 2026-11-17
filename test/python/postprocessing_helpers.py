# Copyright (c) 2026, The Simons Foundation
# This file is part of TRIQS/ctseg and is licensed under the terms of GPLv3 or later.
# SPDX-License-Identifier: GPL-3.0-or-later
# See LICENSE in the root of this distribution for details.

"""Unit tests for the pure helpers in triqs_ctseg.postprocessing.

No Monte Carlo; these tests pin down data-transformation contracts that
are otherwise only exercised indirectly via solve_generic.
"""

from types import SimpleNamespace

import numpy as np
import triqs.utility.mpi as mpi
from triqs.operators import n

from triqs_ctseg.postprocessing import (
    build_color_tables,
    _dd_to_4idx,
    extract_u_tensor_from_h_int,
    check_spectrum,
    compute_sigma_first_moment,
)


def test_build_color_tables_single_orbital():
    # Matches the convention used by anderson.py / dynamical_U.py
    gf_struct = [('down', 1), ('up', 1)]
    blk, idx, orb = build_color_tables(gf_struct)
    assert blk == ['down', 'up']
    assert idx == [0, 0]
    assert orb == [0, 0]


def test_build_color_tables_two_orbital_single_block_per_spin():
    gf_struct = [('up', 2), ('down', 2)]
    blk, idx, orb = build_color_tables(gf_struct)
    assert blk == ['up', 'up', 'down', 'down']
    assert idx == [0, 1, 0, 1]
    assert orb == [0, 1, 0, 1]


def test_build_color_tables_two_orbital_multi_block():
    # One block per (spin, orbital) pair; orbital counter must advance
    # across blocks, independently per spin
    gf_struct = [('up_0', 1), ('up_1', 1), ('down_0', 1), ('down_1', 1)]
    blk, idx, orb = build_color_tables(gf_struct)
    assert blk == ['up_0', 'up_1', 'down_0', 'down_1']
    assert idx == [0, 0, 0, 0]
    assert orb == [0, 1, 0, 1]


def test_dd_to_4idx_2x2():
    direct = np.array([[10., 2.], [2., 20.]])
    exchange = np.array([[0., 0.5], [0.5, 0.]])
    t = _dd_to_4idx(direct, exchange)

    assert t.shape == (2, 2, 2, 2)
    # (i,j,i,j) always set from direct
    assert t[0, 0, 0, 0] == 10.
    assert t[1, 1, 1, 1] == 20.
    assert t[0, 1, 0, 1] == 2.
    assert t[1, 0, 1, 0] == 2.
    # (i,j,j,i) and (i,i,j,j) set from exchange for i != j
    assert t[0, 1, 1, 0] == 0.5
    assert t[1, 0, 0, 1] == 0.5
    assert t[0, 0, 1, 1] == 0.5
    assert t[1, 1, 0, 0] == 0.5
    # Fully-diagonal i==j gets only the direct value (exchange ignored)
    assert t[0, 0, 0, 0] == direct[0, 0]

    # All other entries must remain zero
    for i, j, k, l in np.ndindex(2, 2, 2, 2):
        populated = (
            (i == k and j == l)                    # (i, j, i, j) from direct
            or (i != j and i == l and j == k)      # (i, j, j, i) from exchange
            or (i != k and i == j and k == l)      # (i, i, k, k) from exchange
        )
        if not populated:
            assert t[i, j, k, l] == 0., f"unexpected nonzero at {(i,j,k,l)}: {t[i,j,k,l]}"


def test_dd_to_4idx_3x3():
    n_orb = 3
    direct = np.arange(n_orb * n_orb, dtype=float).reshape(n_orb, n_orb) + 1.
    exchange = np.arange(n_orb * n_orb, dtype=float).reshape(n_orb, n_orb) * 0.1
    t = _dd_to_4idx(direct, exchange)

    for i, j in np.ndindex(n_orb, n_orb):
        assert t[i, j, i, j] == direct[i, j]
        if i != j:
            assert t[i, j, j, i] == exchange[i, j]
            assert t[i, i, j, j] == exchange[i, j]


def test_extract_u_tensor_density_density_two_orbital():
    # Kanamori-like density-density: U (same-orb, opposite-spin),
    # Up (inter-orb opposite-spin), Up - J (inter-orb same-spin)
    U, Up, J = 4.0, 2.0, 0.5
    gf_struct = [('up', 2), ('down', 2)]
    h_int = (
        U * (n('up', 0) * n('down', 0) + n('up', 1) * n('down', 1))
        + Up * (n('up', 0) * n('down', 1) + n('up', 1) * n('down', 0))
        + (Up - J) * (n('up', 0) * n('up', 1) + n('down', 0) * n('down', 1))
    )

    U_dd = extract_u_tensor_from_h_int(h_int, gf_struct, return_4idx=False)
    # Block layout is [up0, up1, down0, down1]
    expected = np.array([
        [0.,        Up - J,    U,         Up],
        [Up - J,    0.,        Up,        U],
        [U,         Up,        0.,        Up - J],
        [Up,        U,         Up - J,    0.],
    ])
    np.testing.assert_allclose(U_dd, expected)

    Uijkl = extract_u_tensor_from_h_int(h_int, gf_struct, return_4idx=True)
    assert Uijkl.shape == (2, 2, 2, 2)
    # Direct = opposite-spin block = [[U, Up], [Up, U]]
    np.testing.assert_allclose(Uijkl[0, 0, 0, 0], U)
    np.testing.assert_allclose(Uijkl[1, 1, 1, 1], U)
    np.testing.assert_allclose(Uijkl[0, 1, 0, 1], Up)
    np.testing.assert_allclose(Uijkl[1, 0, 1, 0], Up)
    # Exchange = direct - same-spin = [[U, Up-(Up-J)], [..., U]] = [[U, J], [J, U]]
    np.testing.assert_allclose(Uijkl[0, 1, 1, 0], J)
    np.testing.assert_allclose(Uijkl[1, 0, 0, 1], J)
    np.testing.assert_allclose(Uijkl[0, 0, 1, 1], J)
    np.testing.assert_allclose(Uijkl[1, 1, 0, 0], J)


def test_check_spectrum_all_within_radius():
    # All eigenvalues |lambda| < 1 -> returned unchanged regardless of truncation flag
    A = np.diag([0.1, 0.5, -0.3])
    np.testing.assert_allclose(check_spectrum(A, radius=1.0, truncation=False), A)
    np.testing.assert_allclose(check_spectrum(A, radius=1.0, truncation=True), A)


def test_check_spectrum_no_truncation_keeps_matrix():
    # Eigenvalue 2.0 is outside radius 1.0 but truncation=False -> returned unchanged
    A = np.diag([0.3, 2.0, -0.1])
    result = check_spectrum(A, radius=1.0, truncation=False)
    np.testing.assert_allclose(result, A)


def test_check_spectrum_truncation_drops_large_eigenvalues():
    # Diagonal so eigendecomposition is trivial; 2.0 and -3.0 should be dropped
    A = np.diag([0.3, 2.0, -0.1, -3.0])
    result = check_spectrum(A, radius=1.0, truncation=True)
    # Surviving eigenvalues must all satisfy |lambda| < 1
    w = np.linalg.eigvals(result)
    assert np.all(np.abs(w) < 1.0 + 1e-10)
    # The two kept eigenvalues (0.3, -0.1) should still be in the spectrum
    w_sorted = np.sort(np.real(w))
    assert np.isclose(w_sorted[0], -0.1) or np.isclose(w_sorted[-1], 0.3)


def test_check_spectrum_rejects_non_square():
    try:
        check_spectrum(np.zeros((2, 3)))
    except ValueError:
        return
    raise AssertionError("check_spectrum should have raised for non-square input")


def _fock_annihilators(n_color):
    """Jordan-Wigner annihilators on the Fock space with state index sum_c n_c 2^c."""
    dim = 2**n_color
    cs = []
    for a in range(n_color):
        c = np.zeros((dim, dim))
        for state in range(dim):
            if state >> a & 1:
                sign = (-1)**bin(state & ((1 << a) - 1)).count('1')
                c[state ^ (1 << a), state] = sign
        cs.append(c)
    return cs


def test_compute_sigma_first_moment_vs_operator_algebra():
    # Kanamori density-density, 2 orbitals, one 2x2 block per spin.
    # Reference: Sigma_0 = -<{[H,c],c+}>, Sigma_1 = <{[H,[H,c]],c+}> - Sigma_0^2
    # evaluated in Fock space for a random diagonal density matrix.
    U, Up, J = 4.0, 2.0, 0.5
    gf_struct = [('up', 2), ('down', 2)]
    h_int = (
        U * (n('up', 0) * n('down', 0) + n('up', 1) * n('down', 1))
        + Up * (n('up', 0) * n('down', 1) + n('up', 1) * n('down', 0))
        + (Up - J) * (n('up', 0) * n('up', 1) + n('down', 0) * n('down', 1))
    )
    # colors [up0, up1, down0, down1]
    U_dd = np.array([
        [0.,     Up - J, U,      Up],
        [Up - J, 0.,     Up,     U],
        [U,      Up,     0.,     Up - J],
        [Up,     U,      Up - J, 0.],
    ])
    n_color = 4

    rng = np.random.default_rng(42)
    prob = rng.random(2**n_color)
    prob /= prob.sum()

    cs = _fock_annihilators(n_color)
    ns = [c.T @ c for c in cs]
    H = 0.5 * sum(U_dd[a, b] * ns[a] @ ns[b] for a in range(n_color) for b in range(n_color))
    rho = np.diag(prob)

    def comm(A, B): return A @ B - B @ A
    def anticomm(A, B): return A @ B + B @ A

    sigma0_ref = np.array([-np.trace(rho @ anticomm(comm(H, c), c.T)) for c in cs])
    sigma1_ref = np.array([np.trace(rho @ anticomm(comm(H, comm(H, c)), c.T)) for c in cs]) - sigma0_ref**2

    # minimal stand-in for the solver; zero D0(tau) / Jperp(tau) mean a static interaction
    solver = SimpleNamespace(
        gf_struct=gf_struct, h_int=h_int,
        results=SimpleNamespace(state_hist=prob),
        D0_tau=[(('up', 'up'), SimpleNamespace(data=np.zeros((3, 2, 2))))],
        Jperp_tau=SimpleNamespace(data=np.zeros((3, 1, 1))),
    )

    sigma1 = compute_sigma_first_moment(solver)
    assert list(sigma1) == ['up', 'down']
    for blk in sigma1.values():
        assert blk.shape == (2, 2)
        np.testing.assert_allclose(blk - np.diag(np.diag(blk)), 0.)
    np.testing.assert_allclose(np.concatenate([np.diag(sigma1['up']), np.diag(sigma1['down'])]),
                               sigma1_ref, rtol=1e-12)

    # no state histogram -> moment unknown
    solver.results.state_hist = None
    assert compute_sigma_first_moment(solver) is None

    # retarded interaction -> moment unknown
    solver.results.state_hist = prob
    solver.Jperp_tau.data[1] = 0.1
    assert compute_sigma_first_moment(solver) is None


if mpi.is_master_node():
    test_build_color_tables_single_orbital()
    test_build_color_tables_two_orbital_single_block_per_spin()
    test_build_color_tables_two_orbital_multi_block()
    test_dd_to_4idx_2x2()
    test_dd_to_4idx_3x3()
    test_extract_u_tensor_density_density_two_orbital()
    test_check_spectrum_all_within_radius()
    test_check_spectrum_no_truncation_keeps_matrix()
    test_check_spectrum_truncation_drops_large_eigenvalues()
    test_check_spectrum_rejects_non_square()
    test_compute_sigma_first_moment_vs_operator_algebra()
    print("postprocessing_helpers: all tests passed")
