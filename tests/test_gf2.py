
import pytest
from homological_compiler.gf2 import GF2Matrix


def test_kernel_and_image():
    A = GF2Matrix([
        [1, 1, 0],
        [0, 1, 1]
    ])

    assert A.rank() == 2
    assert A.kernel_basis() == GF2Matrix([[1], [1], [1]])
    assert (A @ A.kernel_basis()).is_zero()
    assert A.image_basis().rank() == 2


def test_linear_system():
    A = GF2Matrix([
        [1, 1, 0],
        [0, 1, 1]
    ])

    b = GF2Matrix([[1], [0]])
    x = A.solve(b)

    assert A @ x == b


def test_column_space():
    A = GF2Matrix([[1], [0]])

    assert A.in_column_space([1, 0])
    assert not A.in_column_space([0, 1])

    with pytest.raises(ValueError):
        A.solve([0, 1])


def test_empty_matrix():
    A = GF2Matrix.zeros(0, 3)

    assert A.rank() == 0
    assert A.kernel_basis().shape == (3, 3)
