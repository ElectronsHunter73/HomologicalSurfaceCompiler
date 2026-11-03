import itertools

import numpy as np

from homological_compiler.distance import distance
from homological_compiler.layout import SurfaceLayout
from homological_compiler.rotated_surface import Patch


def brute_force(layout, pauli, max_weight):
    code = layout.code()
    checks = (code.hx if pauli == "Z" else code.hz)._data
    x, z = code.logical_operators()
    conjugate = (x if pauli == "Z" else z)._data
    for weight in range(1, max_weight + 1):
        for support in itertools.combinations(range(code.n), weight):
            v = np.zeros(code.n, dtype=np.uint8)
            v[list(support)] = 1
            if not (checks @ v % 2).any() and (conjugate.T @ v % 2).any():
                return weight


def small_layouts():
    two = SurfaceLayout()
    two.add_patch("a", 3)
    two.add_patch("b", 3, next_to="a", side="right")
    stacked = SurfaceLayout()
    stacked.add_patch("a", 3)
    stacked.add_patch("b", 3, next_to="a", side="bottom")
    return [SurfaceLayout([Patch("a", 0, 0, 3, 3)]),
            SurfaceLayout([Patch("a", 0, 0, 5, 5)]),
            SurfaceLayout([Patch("a", 0, 0, 3, 5)]),
            two,
            two.merge("a", "b").after,
            stacked.merge("a", "b").after]


def test_shortest_path_matches_brute_force():
    for layout in small_layouts():
        for pauli in ("Z", "X"):
            assert distance(layout, pauli) == brute_force(layout, pauli, 7)


def test_rectangles():
    for rows, cols in ((3, 7), (5, 9), (7, 3)):
        layout = SurfaceLayout([Patch("a", 0, 0, rows, cols)])
        assert distance(layout, "Z") == rows 
        assert distance(layout, "X") == cols


def test_distance_at_scale():
    d = 15
    layout = SurfaceLayout()
    layout.add_patch("a", d)
    layout.add_patch("b", d, next_to="a", side="right")
    merged = layout.merge("a", "b").after
    assert distance(merged, "Z") == d
    assert distance(merged, "X") == 2 * d + 1
