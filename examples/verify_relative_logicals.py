
import numpy as np

from homological_compiler.chain_complex import rotated_surface_cell_complex
from homological_compiler.rotated_surface import (
    rotated_surface_code,
    qubit_index,
)
from homological_compiler.gf2 import GF2Matrix


d = 5

code = rotated_surface_code(d)
cells = rotated_surface_cell_complex(d)

d1, d2 = cells.relative_boundary_maps()

# Construct a vertical logical Z string.
z_data = np.zeros((code.n, 1), dtype=np.uint8)

for row in range(d):
    qubit = qubit_index(row, d // 2, d)
    z_data[qubit, 0] = 1

z = GF2Matrix(z_data)

# Construct a horizontal logical X string.
x_data = np.zeros((code.n, 1), dtype=np.uint8)

for col in range(d):
    qubit = qubit_index(d // 2, col, d)
    x_data[qubit, 0] = 1

x = GF2Matrix(x_data)

# Verify that Z is a nontrivial relative homology class.
assert (d1 @ z).is_zero()
assert not d2.in_column_space(z)

# Verify that X is a nontrivial relative cohomology class.
assert (d2.T @ x).is_zero()
assert not d1.T.in_column_space(x)

# Verify that the operators anticommute.
assert x.T @ z == GF2Matrix.identity(1)

# Verify the geometric endpoints of Z.
boundary = cells.boundary_1 @ z

assert boundary._data[:-2].sum() == 0
assert boundary._data[-2:, 0].tolist() == [1, 1]

# Compare with the logical basis from Stage 1.
x_basis, z_basis = code.logical_operators()

# Equivalent logical operators differ by stabilizers.
assert code.is_x_stabilizer(x + x_basis)
assert code.is_z_stabilizer(z + z_basis)

print("RELATIVE LOGICAL OPERATORS")
print("Distance:", d)

print("\nLogical Z support:")
print(z._data.reshape(d, d))

print("\nLogical X support:")
print(x._data.reshape(d, d))

print("\nRelative H1 dimension:", cells.relative_h1_dimension())
print("Logical pairing:", (x.T @ z)._data)

print("\nAll logical operator checks passed!")
