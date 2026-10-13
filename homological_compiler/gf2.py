import numpy as np
from sympy import GF
from sympy.polys.matrices import DomainMatrix


_FIELD = GF(2)
class GF2Matrix:

    def __init__(self, data):
        self._data = (np.asarray(data)% 2).astype(np.uint8)

    @classmethod
    def zeros(cls, rows, cols):
        return cls(np.zeros((rows, cols), dtype=np.uint8))

    @classmethod
    def identity(cls, n):
        return cls(np.eye(n, dtype=np.uint8))

    @property
    def shape(self):
        return self._data.shape

    @property
    def n_rows(self):
        return self.shape[0]

    @property
    def n_cols(self):
        return self.shape[1]

    @property
    def T(self):
        return GF2Matrix(self._data.T)

    def is_zero(self):
        return not np.any(self._data)

    def __add__(self, other):
        if not isinstance(other, GF2Matrix):
            return NotImplemented

        if self.shape != other.shape:
            raise ValueError("Matrices must have the same shape")

        return GF2Matrix(self._data ^ other._data)

    def __sub__(self, other):
        return self.__add__(other)

    def __matmul__(self, other):
        if not isinstance(other, GF2Matrix):
            return NotImplemented

        if self.n_cols != other.n_rows:
            raise ValueError("Incompatible shapes for matrix multiplication")

        return GF2Matrix((self._data @ other._data) % 2)

    def __eq__(self, other):
        if not isinstance(other, GF2Matrix):
            return NotImplemented

        return np.array_equal(self._data, other._data)

    def __repr__(self):
        return f"GF2Matrix({self._data.tolist()!r})"

    def _domain_matrix(self):
        if self.n_rows == 0:
            return DomainMatrix.zeros(self.shape, _FIELD)

        return DomainMatrix.from_list(
            self._data.tolist(), _FIELD
        )

    @staticmethod
    def _from_domain_matrix(matrix):
        data = np.array(
            [
                [int(value) for value in row]
                for row in matrix.to_list()
            ],
            dtype=np.uint8,
        ).reshape(matrix.shape)

        return GF2Matrix(data)

    def rref(self):
        reduced, pivots = self._domain_matrix().rref()

        return self._from_domain_matrix(reduced), tuple(pivots)

    def rank(self):
        return len(self.rref()[1])

    def kernel_basis(self):
        basis_rows = self._domain_matrix().nullspace()

        return self._from_domain_matrix(basis_rows).T

    def image_basis(self):
        pivots = self.rref()[1]

        return GF2Matrix(self._data[:, list(pivots)])

    def solve(self, rhs):
        if isinstance(rhs, GF2Matrix):
            b = rhs
        else:
            data = np.asarray(rhs)
            b = GF2Matrix(
                data[:, None] if data.ndim == 1 else data
            )

        if b.n_rows != self.n_rows:
            raise ValueError("Wrong right-hand side dimension")

        combined = GF2Matrix(
            np.column_stack((self._data, b._data))
        )

        reduced, pivots = combined.rref()

        if any(pivot >= self.n_cols for pivot in pivots):
            raise ValueError("The linear system has no solution")

        x = np.zeros(
            (self.n_cols, b.n_cols),
            dtype=np.uint8
        )

        for row, col in enumerate(pivots):
            x[col] = reduced._data[row, self.n_cols:]

        return GF2Matrix(x)

    def in_column_space(self, vector):
        if isinstance(vector, GF2Matrix):
            v = vector
        else:
            data = np.asarray(vector)
            v = GF2Matrix(
                data[:, None] if data.ndim == 1 else data
            )

        if v.shape != (self.n_rows, 1):
            raise ValueError(
                "Expected a single column vector of matching length"
            )

        combined = GF2Matrix(
            np.column_stack((self._data, v._data))
        )

        return self.rank() == combined.rank()
