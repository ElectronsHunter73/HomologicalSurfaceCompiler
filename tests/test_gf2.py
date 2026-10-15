import pytest
from homological_compiler.gf2 import GF2Matrix
import numpy as np

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
def test_rank_nullity():
    rng = np.random.default_rng(42)

    for m in range(5):
        for n in range(5):
            A = GF2Matrix(rng.integers(0, 2, size=(m, n)))

            K = A.kernel_basis()

            assert A.rank() + K.n_cols == n
            assert K.rank() == K.n_cols
            assert (A @ K).is_zero()
            assert A.image_basis().rank() == A.rank()


def test_matrix_multiplication():
    A = GF2Matrix([[1, 1], [0, 1]])
    B = GF2Matrix([[1, 0], [1, 1]])

    assert A @ B == GF2Matrix([[0, 1], [1, 1]])

    A = GF2Matrix(np.ones((1, 513), dtype=np.uint8))
    B = GF2Matrix(np.ones((513, 1), dtype=np.uint8))

    assert A @ B == GF2Matrix([[1]])