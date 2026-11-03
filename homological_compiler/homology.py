import numpy as np
from .gf2 import GF2Matrix


class Prediction:

    def __init__(self, surgery):
        self.surgery = surgery
        if surgery.kind == "merge":
            self.small, self.big = surgery.before, surgery.after
        else:
            self.small, self.big = surgery.after, surgery.before

        Y, X = self.small.complex(), self.big.complex()
        self.pauli = "Z"
        if not Y.is_subcomplex_of(X):
            Y, X = Y.dual(), X.dual()
            self.pauli = "X"
            if not Y.is_subcomplex_of(X):
                raise ValueError("Neither complex of the smaller layout includes into the larger")
        conjugate = "X" if self.pauli == "Z" else "Z"
        self.relative = Q = X.relative_to(Y)
        self.h2 = Q.cycles(2)
        self.generators = [Q.support(2, self.h2._data[:, j]) for j in range(self.h2.n_cols)]
        self.measured = [X.boundary(2, S) for S in self.generators]

        patches = self.small.patches
        pairing = np.zeros((len(self.measured), len(patches)), dtype=np.uint8)
        for i, op in enumerate(self.measured):
            for j, p in enumerate(patches):
                pairing[i, j] = len(op & p.logical(conjugate)) % 2
        self.labels = [" ".join(f"{self.pauli}[{p.name}]" for p, bit in zip(patches, row) if bit)
                       for row in pairing]

        names = surgery.inputs + surgery.outputs
        involved = [j for j, p in enumerate(patches) if p.name in names]
        K = GF2Matrix(pairing[:, involved]).kernel_basis()
        self.survivors = [" ".join(f"{conjugate}[{patches[involved[i]].name}]"
                                   for i in np.flatnonzero(K._data[:, j]))
                          for j in range(K.n_cols)]

        k_small = len(patches)
        k_big = k_small - GF2Matrix(pairing).rank() + Q.betti(1)
        if surgery.kind == "merge":
            self.k_before, self.k_after = k_small, k_big
        else:
            self.k_before, self.k_after = k_big, k_small

    def report(self):
        s, Q = self.surgery, self.relative
        complex_name = "primal" if self.pauli == "Z" else "dual"
        lines = [
            f"{s.kind} {' + '.join(s.inputs)} -> {', '.join(s.outputs)}: "
            f"inclusion of {complex_name} complexes",
            f"  C(X,Y): {len(Q.cells[1])} strip qubits, {len(Q.cells[2])} strip checks, "
            f"{len(Q.cells[0])} new vertices",
        ]
        for label, S in zip(self.labels, self.generators):
            if s.kind == "merge":
                lines.append(f"  measures {label} = product of {len(S)} new checks")
            else:
                lines.append(f"  lifts differ by {label}, whose value the merge fixed")
        if s.kind == "merge":
            lines.append(f"  survives: {', '.join(self.survivors)}")
        lines.append(f"  logical qubits: {self.k_before} -> {self.k_after}")
        return "\n".join(lines)


def gauge_fix(prediction, outcomes):
    Q, S = prediction.relative, prediction.h2
    f = Q.chain(2, [site for site, bit in outcomes.items() if bit])
    m = S.T @ f
    E = S.T.solve(GF2Matrix.identity(S.n_cols))
    c = Q.d[2].T.solve(f + E @ m)
    return [int(bit) for bit in m._data[:, 0]], Q.support(1, c._data)
