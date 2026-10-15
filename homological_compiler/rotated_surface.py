import numpy as np
from .css_code import CSSCode
from .gf2 import GF2Matrix


def qubit_index(row, col, d):
    if (
        not isinstance(row, int)
        or not isinstance(col, int)
        or not (0 <= row < d and 0 <= col < d)
    ):
        raise ValueError("Invalid qubit coordinates")

    return row * d + col


def interior_checks(d):
    x_checks = []
    z_checks = []

    for row in range(d - 1):
        for col in range(d - 1):
            support = [
                qubit_index(row, col, d),
                qubit_index(row, col + 1, d),
                qubit_index(row + 1, col, d),
                qubit_index(row + 1, col + 1, d),
            ]

            if (row + col) % 2 == 0:
                x_checks.append(support)
            else:
                z_checks.append(support)

    return x_checks, z_checks


def boundary_checks(d):
    x_checks = []
    z_checks = []

    for row in range(1, d - 1, 2):
        support = [qubit_index(row, 0, d),qubit_index(row + 1, 0, d)]
        x_checks.append(support)
    for row in range(0, d - 1, 2):
        support = [qubit_index(row, d - 1, d),qubit_index(row + 1, d - 1, d)]
        x_checks.append(support)

    for col in range(0, d - 1, 2):
        support = [qubit_index(0, col, d),qubit_index(0, col + 1, d)]
        z_checks.append(support)

    for col in range(1, d - 1, 2):
        support = [qubit_index(d - 1, col, d),qubit_index(d - 1, col + 1, d)]
        z_checks.append(support)

    return x_checks, z_checks


def _checks_to_matrix(checks, n):
    data = np.zeros((len(checks), n),dtype=np.uint8)
    for i, support in enumerate(checks):
        data[i, support] = 1
    return GF2Matrix(data)


def rotated_surface_code(d):
    interior_x, interior_z = interior_checks(d)
    boundary_x, boundary_z = boundary_checks(d)
    x_checks = interior_x + boundary_x
    z_checks = interior_z + boundary_z
    n = d * d
    hx = _checks_to_matrix(x_checks, n)
    hz = _checks_to_matrix(z_checks, n)
    return CSSCode(hx,hz)

