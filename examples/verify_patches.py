
from homological_compiler.patch_collection import PatchCollection
from homological_compiler.gf2 import GF2Matrix


distances = [3, 5, 7]

patches = PatchCollection()

for d in distances:
    patches.add_patch(d)

cells = patches.cell_complex()
code = patches.css_code()

print("MULTIPLE SURFACE-CODE PATCHES")
print("Distances:", distances)
print("Number of patches:", len(patches.patches))

print("\nCOMBINED CELL COMPLEX")
print("Vertices:", len(cells.vertices))
print("Edges:", len(cells.edges))
print("Faces:", len(cells.faces))
print("Boundary vertices:", cells.boundary_vertices)

# Verify the chain-complex condition.
assert (cells.boundary_1 @ cells.boundary_2).is_zero()

# Verify the number of physical and logical qubits.
assert code.n == sum(d * d for d in distances)
assert code.k == len(distances)

# Verify relative homology.
assert cells.relative_h1_dimension() == code.k

print("\nCOMBINED CSS CODE")
print("Physical qubits:", code.n)
print("Logical qubits:", code.k)
print("Relative H1 dimension:", cells.relative_h1_dimension())

# Verify the logical operators.
x, z = code.logical_operators()

assert (code.hz @ x).is_zero()
assert (code.hx @ z).is_zero()

assert x.T @ z == GF2Matrix.identity(code.k)

print("\nLogical pairing:")
print((x.T @ z)._data)

print("\nAll multiple-patch checks passed!")
