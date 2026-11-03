# HomologicalSurfaceCompiler

**Work in progress.** This is not a finished compiler yet. The parts that work
are listed below, and the rest is still on the todo list.

Lattice surgery on rotated surface codes in Python. You place surface code
patches on a grid and merge and split them. For every step the code works out
what happens to the logical qubits: which logical operator is measured, how to
read its outcome, which logical operators are still valid afterwards, and which
correction to apply. It gets this from the homology of the strip between the
patches, and checks every answer against a direct GF(2) computation on the
stabilizer matrices.

## What works

- surface code patches and layouts on a grid
- ZZ and XX merges and splits between neighbouring patches
- for each merge: the measured logical operator, its outcome from the new
  checks, the logical operators that survive, and the correction on the strip
- the number of logical qubits after each step
- GF(2) verification of all of the above
- code distance (as a shortest path)

## Not done yet

- measuring a single patch
- following the logical operators through a whole program (Pauli frame), for
  example combining the outcomes below into the corrections of a CNOT
- export to Stim and OpenQASM, so no circuit-level simulation or noise yet
- qLDPC codes

## Install

```
git clone https://github.com/ElectronsHunter73/HomologicalSurfaceCompiler.git
cd HomologicalSurfaceCompiler
pip install -e ".[dev]"
```

Needs Python 3.10+, numpy and sympy.

## Example: the measurements of a lattice-surgery CNOT

A CNOT by lattice surgery uses three patches: a control C, a target T and an
ancilla A prepared in |+>. You measure Z[C] Z[A], then X[A] X[T], then A on its
own. The code below does the first two measurements as merges and splits. C
sits left of A, so their merge is a ZZ merge. T sits below A, so their merge is
an XX merge. The code only gets the layout and works out the rest. The new
checks give random results, like on hardware, so here they are drawn at random.

```python
import random

from homological_compiler.layout import SurfaceLayout
from homological_compiler.homology import Prediction, gauge_fix
from homological_compiler.verification import verify

random.seed(1)

layout = SurfaceLayout()
layout.add_patch("A", 3)
layout.add_patch("C", 3, next_to="A", side="left")
layout.add_patch("T", 3, next_to="A", side="bottom")
print(layout.draw())

def step(surgery):
    prediction = Prediction(surgery)
    assert verify(prediction) == []
    print()
    print(prediction.report())
    if surgery.kind == "merge":
        outcomes = {site: random.randint(0, 1) for site in surgery.added}
        m, correction = gauge_fix(prediction, outcomes)
        print("  new check outcomes:", list(outcomes.values()))
        print("  logical outcome m =", m[0], "| correction on the strip:", sorted(correction))
    return surgery.after

layout = step(layout.merge("C", "A"))
layout = step(layout.split("C+A"))
layout = step(layout.merge("A", "T"))
layout = step(layout.split("A+T"))
```

Output:

```
    Z               Z
  o   o   o       o   o   o
    X   Z   X       X   Z   X
  o   o   o       o   o   o
X   Z   X       X   Z   X
  o   o   o       o   o   o
        Z               Z

                    Z
                  o   o   o
                    X   Z   X
                  o   o   o
                X   Z   X
                  o   o   o
                        Z

merge C + A -> C+A: inclusion of primal complexes
  C(X,Y): 3 strip qubits, 4 strip checks, 0 new vertices
  measures Z[A] Z[C] = product of 4 new checks
  survives: X[A] X[C]
  logical qubits: 3 -> 2
  new check outcomes: [0, 0, 1, 0]
  logical outcome m = 1 | correction on the strip: [q(0,-1), q(1,-1)]

split C+A -> C, A: inclusion of primal complexes
  C(X,Y): 3 strip qubits, 4 strip checks, 0 new vertices
  lifts differ by Z[C] Z[A], whose value the merge fixed
  logical qubits: 2 -> 3

merge A + T -> A+T: inclusion of dual complexes
  C(X,Y): 3 strip qubits, 4 strip checks, 0 new vertices
  measures X[T] X[A] = product of 4 new checks
  survives: Z[T] Z[A]
  logical qubits: 3 -> 2
  new check outcomes: [1, 1, 1, 1]
  logical outcome m = 0 | correction on the strip: [q(3,0), q(3,2)]

split A+T -> A, T: inclusion of dual complexes
  C(X,Y): 3 strip qubits, 4 strip checks, 0 new vertices
  lifts differ by X[A] X[T], whose value the merge fixed
  logical qubits: 2 -> 3
```

In the first merge the four new checks give 0, 0, 1, 0, so the measured value
of Z[A] Z[C] is m = 1. X[A] X[C] is still a logical operator of the merged
patch, and X on the two listed strip qubits brings the merged patch back to a
clean code state. The split then gives back C and A.

There is also a demo with four patches, drawings of each step and the code
distance:

```
python examples/surgery_demo.py
```



## Files

```
homological_compiler/
    gf2.py              linear algebra over GF(2)
    css_code.py         CSS codes
    geometry.py         qubits and check sites on the grid
    rotated_surface.py  one surface code patch
    layout.py           layouts, merge and split
    chain_complex.py    chain complexes over GF(2)
    homology.py         predictions from homology
    verification.py     GF(2) checks of the predictions
    distance.py         code distance
examples/
tests/
```
