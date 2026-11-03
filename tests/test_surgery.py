import numpy as np

from homological_compiler.homology import Prediction, gauge_fix
from homological_compiler.layout import SurfaceLayout
from homological_compiler.verification import verify


def grid(d):
    layout = SurfaceLayout()
    layout.add_patch("p0", d)
    layout.add_patch("p1", d, next_to="p0", side="right")
    layout.add_patch("p2", d, next_to="p0", side="bottom")
    layout.add_patch("p3", d, next_to="p1", side="bottom")
    return layout


def test_merge_predictions_are_confirmed_by_gf2():
    for d in (3, 5, 7):
        for a, b, measured, survivor in (("p0", "p1", "Z[p0] Z[p1]", "X[p0] X[p1]"),
                                         ("p0", "p2", "X[p0] X[p2]", "Z[p0] Z[p2]")):
            prediction = Prediction(grid(d).merge(a, b))
            assert prediction.labels == [measured]
            assert prediction.survivors == [survivor]
            assert (prediction.k_before, prediction.k_after) == (4, 3)
            assert len(prediction.generators[0]) == d + 1
            assert verify(prediction) == []


def test_split_predictions_are_confirmed_by_gf2():
    for d in (3, 5):
        merged = grid(d).merge("p0", "p1").after
        prediction = Prediction(merged.split("p0+p1"))
        assert prediction.labels == ["Z[p0] Z[p1]"]
        assert (prediction.k_before, prediction.k_after) == (3, 4)
        assert verify(prediction) == []


def test_a_program_of_merges_and_splits():
    layout = grid(3)
    program = [("merge", "p0", "p1", "Z[p0] Z[p1]"),
               ("merge", "p2", "p3", "Z[p2] Z[p3]"),
               ("merge", "p0+p1", "p2+p3", "X[p0+p1] X[p2+p3]"),
               ("split", "p0+p1+p2+p3", None, "X[p0+p1] X[p2+p3]"),
               ("split", "p0+p1", None, "Z[p0] Z[p1]"),
               ("split", "p2+p3", None, "Z[p2] Z[p3]")]
    current = layout
    for kind, a, b, label in program:
        surgery = current.merge(a, b) if kind == "merge" else current.split(a)
        prediction = Prediction(surgery)
        assert prediction.labels == [label]
        assert verify(prediction) == []
        current = surgery.after
    assert current.checks() == layout.checks()


def test_a_row_of_three_patches():
    layout = SurfaceLayout()
    layout.add_patch("a", 3)
    layout.add_patch("b", 3, next_to="a", side="right")
    layout.add_patch("c", 3, next_to="b", side="right")
    prediction = Prediction(layout.merge("a", "b").after.merge("a+b", "c"))
    assert prediction.labels == ["Z[a+b] Z[c]"]
    assert verify(prediction) == []


def test_gauge_fix_separates_the_logical_result_from_gauge():
    rng = np.random.default_rng(7)
    for a, b in (("p0", "p1"), ("p0", "p2")):
        prediction = Prediction(grid(5).merge(a, b))
        Q, S = prediction.relative, prediction.generators[0]
        for _ in range(20):
            outcomes = {site: int(rng.integers(2)) for site in Q.cells[2]}
            m, correction = gauge_fix(prediction, outcomes)
            assert m == [sum(outcomes[s] for s in S) % 2]
            assert correction <= prediction.surgery.strip
            flipped = Q.support(2, (Q.d[2].T @ Q.chain(1, correction))._data)
            remaining = {site for site, bit in outcomes.items() if bit} ^ flipped
            assert len(remaining) == m[0]
