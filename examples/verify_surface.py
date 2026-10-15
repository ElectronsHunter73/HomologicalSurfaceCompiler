
from homological_compiler.rotated_surface import rotated_surface_code
from homological_compiler.gf2 import GF2Matrix

d = 5

code = rotated_surface_code(d)

print("Physical qubits:", code.n)
print("Logical qubits:", code.k)
print("X stabilizer rank:", code.rank_x)
print("Z stabilizer rank:", code.rank_z)

assert (code.hx @ code.hz.T).is_zero()

assert code.n == d**2
assert code.k == 1
assert code.rank_x + code.rank_z == d**2 - 1

x, z = code.logical_operators()

assert (code.hz @ x).is_zero()
assert (code.hx @ z).is_zero()
assert x.T @ z == GF2Matrix.identity(1)

print("\nLogical X:", x)
print("Logical Z:", z)
print("\nAll algebraic checks passed!")
