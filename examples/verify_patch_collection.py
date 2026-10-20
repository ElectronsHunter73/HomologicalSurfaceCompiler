
import numpy as np

from homological_compiler.patch_collection import PatchCollection
from homological_compiler.gf2 import GF2Matrix


for d in (3, 5):

    patches = PatchCollection()

    p0 = patches.add_patch(d)
    p1 = patches.add_patch(d, next_to=p0, side="right")
    p2 = patches.add_patch(d, next_to=p0, side="bottom")
    p3 = patches.add_patch(d, next_to=p1, side="bottom")

    code = patches.css_code()

    assert code.k == 4

    print(f"\nDISTANCE {d}")

    for a, b, side_a, side_b in patches.adjacent_pairs():

        measurement, pairs = patches.candidate_interface(a, b)

        data = np.zeros((code.n, 1), dtype=np.uint8)

        for qa, qb in pairs:
            data[qa, 0] = 1
            data[qb, 0] = 1

        operator = GF2Matrix(data)

        if measurement == "ZZ":
            assert code.is_z_logical(operator)
            dual_side = "top"

        else:
            assert code.is_x_logical(operator)
            dual_side = "left"

        for i in range(len(patches.patches)):

            dual_data = np.zeros((code.n, 1), dtype=np.uint8)

            for q in patches.boundary_qubits(i, dual_side):
                dual_data[q, 0] = 1

            dual = GF2Matrix(dual_data)

            pairing = (dual.T @ operator)._data[0, 0]

            assert pairing == int(i in (a, b))

        print(f"Interface {a}-{b}: {measurement} verified")
        

    print("All logical interface checks passed!")
