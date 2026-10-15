import numpy as np

from .gf2 import GF2Matrix


class CSSCode:
    def __init__(self, hx, hz):
        if not isinstance(hx, GF2Matrix) or not isinstance(hz, GF2Matrix):
            raise TypeError("hx and hz must be GF2Matrix objects")

        if hx.n_cols != hz.n_cols:
            raise ValueError("Hx and Hz must have the same number of columns")

        if not (hx @ hz.T).is_zero():
            raise ValueError("CSS checks must commute: Hx @ Hz.T != 0")

        self.hx = hx
        self.hz = hz

    @property
    def n(self):
        return self.hx.n_cols

    @property
    def rank_x(self):
        return self.hx.rank()

    @property
    def rank_z(self):
        return self.hz.rank()

    @property
    def k(self):
        return self.n - self.rank_x - self.rank_z

    def _as_column(self, operator):
        if isinstance(operator, GF2Matrix):
            vector = operator
        else:
            data = np.asarray(operator)
            vector = GF2Matrix(data[:, None] if data.ndim == 1 else data)

        if vector.shape != (self.n, 1):
            raise ValueError(
                f"Expected an operator vector of shape ({self.n}, 1)"
            )

        return vector

    def z_syndrome(self, z):
        return self.hx @ self._as_column(z)

    def x_syndrome(self, x):
        return self.hz @ self._as_column(x)

    def is_z_stabilizer(self, z):
        return self.hz.T.in_column_space(self._as_column(z))

    def is_x_stabilizer(self, x):
        return self.hx.T.in_column_space(self._as_column(x))

    def is_z_logical(self, z):
        vector = self._as_column(z)
        return (
            (self.hx @ vector).is_zero()
            and not self.is_z_stabilizer(vector)
        )

    def is_x_logical(self, x):
        vector = self._as_column(x)
        return (
            (self.hz @ vector).is_zero()
            and not self.is_x_stabilizer(vector)
        )

    @staticmethod
    def _quotient_basis(kernel: GF2Matrix,trivial: GF2Matrix,dimension: int):
        span = trivial
        chosen = []

        for j in range(kernel.n_cols):
            candidate = GF2Matrix(kernel._data[:, j:j + 1])

            if not span.in_column_space(candidate):
                chosen.append(candidate._data)

                span = GF2Matrix(
                    np.column_stack((span._data, candidate._data))
                )

                if len(chosen) == dimension:
                    break

        if len(chosen) != dimension:
            raise RuntimeError(
                "Could not construct the expected logical basis"
            )

        if not chosen:
            return GF2Matrix.zeros(kernel.n_rows, 0)

        return GF2Matrix(np.column_stack(chosen))

    def logical_z_basis(self):
        return self._quotient_basis(
            self.hx.kernel_basis(),
            self.hz.T,
            self.k
        )

    def logical_x_basis(self):
        return self._quotient_basis(
            self.hz.kernel_basis(),
            self.hx.T,
            self.k
        )

    def logical_operators(self) -> tuple[GF2Matrix, GF2Matrix]:
        x = self.logical_x_basis()
        z = self.logical_z_basis()

        if self.k == 0:
            return x, z

        pairing = x.T @ z

        inverse = pairing.solve(GF2Matrix.identity(self.k))

        return x, z @ inverse
