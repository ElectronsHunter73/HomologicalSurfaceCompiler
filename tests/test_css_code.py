import pytest

from homological_compiler.css_code import CSSCode
from homological_compiler.gf2 import GF2Matrix


def four_qubit_code():
    # A valid [[4, 1]] CSS code used by several tests.
    hx = GF2Matrix([[1, 1, 1, 1]])

    hz = GF2Matrix([
        [1, 1, 0, 0],
        [0, 0, 1, 1]
    ])

    return CSSCode(hx, hz)


def test_valid_code():
    # Four physical qubits and one logical qubit.
    code = four_qubit_code()

    assert code.n == 4
    assert code.rank_x == 1
    assert code.rank_z == 2
    assert code.k == 1


def test_noncommuting_checks():
    # X and Z overlap on exactly one qubit.
    # They anticommute, so the CSS code is invalid.
    hx = GF2Matrix([[1, 0]])
    hz = GF2Matrix([[1, 0]])

    with pytest.raises(ValueError):
        CSSCode(hx, hz)


def test_different_qubit_counts():
    # Hx and Hz must refer to the same physical qubits.
    hx = GF2Matrix([[1, 1]])
    hz = GF2Matrix([[1, 1, 0]])

    with pytest.raises(ValueError):
        CSSCode(hx, hz)


def test_wrong_input_types():
    # Both inputs must be GF2Matrix objects.
    hx = [[1, 1]]
    hz = GF2Matrix([[1, 1]])

    with pytest.raises(TypeError):
        CSSCode(hx, hz)


def test_dependent_checks_do_not_change_k():
    # Repeating a stabilizer does not add an independent
    # constraint, so the number of logical qubits is unchanged.
    hx = GF2Matrix([
        [1, 1, 1, 1],
        [1, 1, 1, 1]
    ])

    hz = GF2Matrix([
        [1, 1, 0, 0],
        [0, 0, 1, 1]
    ])

    code = CSSCode(hx, hz)

    assert code.rank_x == 1
    assert code.rank_z == 2
    assert code.k == 1


def test_syndromes():
    code = four_qubit_code()

    # A Z error on qubit 0 anticommutes with the X check.
    assert code.z_syndrome([1, 0, 0, 0]) == GF2Matrix([
        [1]
    ])

    # An X error on qubit 0 anticommutes with the first
    # Z check, but commutes with the second.
    assert code.x_syndrome([1, 0, 0, 0]) == GF2Matrix([
        [1],
        [0]
    ])


def test_stabilizers_and_logicals():
    code = four_qubit_code()

    # These vectors describe products of stabilizer generators.
    assert code.is_z_stabilizer([1, 1, 0, 0])
    assert code.is_x_stabilizer([1, 1, 1, 1])

    # These commute with every stabilizer, but they
    # are not themselves products of stabilizers.
    assert code.is_z_logical([1, 0, 1, 0])
    assert code.is_x_logical([1, 1, 0, 0])

    # A stabilizer is not a nontrivial logical operator.
    assert not code.is_z_logical([1, 1, 0, 0])

    # An error with nonzero syndrome is not logical.
    assert not code.is_z_logical([1, 0, 0, 0])
    assert not code.is_x_logical([1, 0, 0, 0])

    # The identity is not a nontrivial logical operator.
    assert not code.is_z_logical([0, 0, 0, 0])


def test_operator_shape_validation():
    # A Pauli operator must have one entry per physical qubit.
    code = four_qubit_code()

    with pytest.raises(ValueError):
        code.z_syndrome([1, 0])


def test_logical_bases():
    code = four_qubit_code()

    x, z = code.logical_operators()

    # One independent logical X and one logical Z.
    assert x.shape == (4, 1)
    assert z.shape == (4, 1)

    # Logical X must commute with all Z stabilizers.
    assert (code.hz @ x).is_zero()

    # Logical Z must commute with all X stabilizers.
    assert (code.hx @ z).is_zero()

    # Logical X and logical Z must anticommute:
    # their binary overlap is 1.
    assert x.T @ z == GF2Matrix.identity(1)

    assert code.is_x_logical(x)
    assert code.is_z_logical(z)


def test_zero_x_checks():
    # A CSS code may have no X stabilizers.
    hx = GF2Matrix.zeros(0, 3)
    hz = GF2Matrix([[1, 1, 0]])

    code = CSSCode(hx, hz)

    assert code.k == 2

    x, z = code.logical_operators()

    assert x.shape == (3, 2)
    assert z.shape == (3, 2)

    # The two logical X/Z pairs must have canonical pairing.
    assert x.T @ z == GF2Matrix.identity(2)


def test_zero_logicals():
    # Both independent stabilizers constrain the two
    # physical qubits, leaving no logical qubits.
    hx = GF2Matrix([[1, 1]])
    hz = GF2Matrix([[1, 1]])

    code = CSSCode(hx, hz)

    assert code.k == 0

    x, z = code.logical_operators()

    # There are no nontrivial logical basis vectors.
    assert x.shape == (2, 0)
    assert z.shape == (2, 0)


def test_two_logical_qubits():
    # Six physical qubits, one independent X check
    # and three independent Z checks.
    #
    # k = 6 - 1 - 3 = 2
    hx = GF2Matrix([
        [1, 1, 1, 1, 0, 0]
    ])

    hz = GF2Matrix([
        [1, 1, 0, 0, 0, 0],
        [0, 0, 1, 1, 0, 0],
        [0, 0, 0, 0, 1, 1]
    ])

    code = CSSCode(hx, hz)

    assert code.k == 2

    # Compute independent representatives of the
    # logical X and Z equivalence classes.
    x = code.logical_x_basis()
    z = code.logical_z_basis()

    assert x.shape == (6, 2)
    assert z.shape == (6, 2)

    # They must preserve the code space.
    assert (hz @ x).is_zero()
    assert (hx @ z).is_zero()

    # Canonical pairing of the two logical qubits.
    paired_x, paired_z = code.logical_operators()

    assert paired_x.T @ paired_z == GF2Matrix.identity(2)


def test_no_stabilizers():
    # Without stabilizers, every physical qubit
    # is also a logical qubit.
    hx = GF2Matrix.zeros(0, 3)
    hz = GF2Matrix.zeros(0, 3)

    code = CSSCode(hx, hz)

    assert code.n == 3
    assert code.rank_x == 0
    assert code.rank_z == 0
    assert code.k == 3

    x, z = code.logical_operators()

    assert x.shape == (3, 3)
    assert z.shape == (3, 3)

    assert x.T @ z == GF2Matrix.identity(3)


def test_stabilizers_have_zero_syndrome():
    # Stabilizers commute with all stabilizers,
    # so their corresponding syndromes must vanish.
    code = four_qubit_code()

    assert code.z_syndrome([1, 1, 0, 0]).is_zero()
    assert code.x_syndrome([1, 1, 1, 1]).is_zero()


def test_reject_wrong_shape_for_other_methods():
    code = four_qubit_code()

    with pytest.raises(ValueError):
        code.is_z_stabilizer([1, 0])

    with pytest.raises(ValueError):
        code.is_x_logical([1, 0, 0])

    # A row matrix is not an n-by-1 column vector.
    with pytest.raises(ValueError):
        code.x_syndrome([[1, 0, 0, 0]])
