
import numpy as np

from homological_compiler.chain_complex import rotated_surface_cell_complex
from homological_compiler.rotated_surface import rotated_surface_code
from homological_compiler.gf2 import GF2Matrix


d = 5

code = rotated_surface_code(d)
cells = rotated_surface_cell_complex(d)

print("GEOMETRIC CELL COMPLEX")
print("Distance:", d)
print("Vertices:", len(cells.vertices))
print("Edges:", len(cells.edges))
print("Faces:", len(cells.faces))
print("Boundary vertices:", cells.boundary_vertices)

assert (cells.boundary_1 @ cells.boundary_2).is_zero()

assert GF2Matrix(cells.boundary_1._data[:-2, :]) == code.hx

assert cells.boundary_2 == code.hz.T

print("\nCHAIN COMPLEX")
print("Boundary 1 shape:", cells.boundary_1.shape)
print("Boundary 2 shape:", cells.boundary_2.shape)

h1 = (
    len(cells.edges)
    - cells.boundary_1.rank()
    - cells.boundary_2.rank()
)

print("Absolute H1 dimension:", h1)

assert h1 == 0

data = np.zeros((code.n, 1), dtype=np.uint8)

for row in range(d):
    qubit = row * d + d // 2
    data[qubit, 0] = 1

z = GF2Matrix(data)

assert code.is_z_logical(z)

assert (code.hx @ z).is_zero()

boundary = cells.boundary_1 @ z

assert boundary._data[:-2].sum() == 0
assert boundary._data[-2:, 0].tolist() == [1, 1]
print("\nLOGICAL Z STRING")
print("Top endpoint:", int(boundary._data[-2, 0]))
print("Bottom endpoint:", int(boundary._data[-1, 0]))
print("\nAll geometric checks passed!")
d1, d2 = cells.relative_boundary_maps()
assert (d1 @ d2).is_zero()
assert d1 == code.hx
assert d2 == code.hz.T
print("\nRelative H1 dimension:", cells.relative_h1_dimension())
assert cells.relative_h1_dimension() == code.k == 1
