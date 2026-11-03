import pytest

from homological_compiler.layout import SurfaceLayout
from homological_compiler.rotated_surface import Patch, rotated_surface_code


def grid(d):

    layout = SurfaceLayout()
    layout.add_patch("p0", d)
    layout.add_patch("p1", d, next_to="p0", side="right")
    layout.add_patch("p2", d, next_to="p0", side="bottom")
    layout.add_patch("p3", d, next_to="p1", side="bottom")
    return layout


def original_checks(d):
    x, z = set(), set()
    for r in range(d - 1):
        for c in range(d - 1):
            square = frozenset({(r, c), (r, c + 1), (r + 1, c), (r + 1, c + 1)})
            (x if (r + c) % 2 == 0 else z).add(square)
    for r in range(1, d - 1, 2):
        x.add(frozenset({(r, 0), (r + 1, 0)}))
    for r in range(0, d - 1, 2):
        x.add(frozenset({(r, d - 1), (r + 1, d - 1)}))
    for c in range(0, d - 1, 2):
        z.add(frozenset({(0, c), (0, c + 1)}))
    for c in range(1, d - 1, 2):
        z.add(frozenset({(d - 1, c), (d - 1, c + 1)}))
    return x, z


def rows_as_sets(matrix, d):
    return {frozenset(divmod(int(j), d) for j in row.nonzero()[0]) for row in matrix._data}


def test_single_patch_is_the_original_code():
    for d in (3, 5, 7, 9):
        code = rotated_surface_code(d)
        x, z = original_checks(d)
        assert rows_as_sets(code.hx, d) == x
        assert rows_as_sets(code.hz, d) == z
        assert code.k == 1


def test_every_rectangle_encodes_one_qubit():
    for rows in range(2, 7):
        for cols in range(2, 7):
            for row in (0, 1): 
                patch = Patch("p", row, 0, rows, cols)
                layout = SurfaceLayout([patch])
                code = layout.code()
                assert code.k == 1
                assert code.rank_x == code.hx.n_rows and code.rank_z == code.hz.n_rows
                assert code.is_z_logical(layout.vector(patch.logical("Z")))
                assert code.is_x_logical(layout.vector(patch.logical("X")))
                assert len(patch.logical("Z") & patch.logical("X")) == 1


def test_patches_without_a_strip_are_rejected():
    layout = SurfaceLayout([Patch("a", 0, 0, 3, 3), Patch("b", 0, 3, 3, 3)])
    with pytest.raises(ValueError):
        layout.checks()


def test_neighbours():
    assert sorted(grid(3).neighbours()) == [
        ("p0", "p1", "ZZ"), ("p0", "p2", "XX"), ("p1", "p3", "XX"), ("p2", "p3", "ZZ"),
    ]


def test_merge_changes_only_the_seam():
    for d in (3, 5, 7):
        for a, b, measurement in (("p0", "p1", "ZZ"), ("p0", "p2", "XX")):
            surgery = grid(d).merge(a, b)
            new_type = measurement[0]
            old = surgery.before.checks()
            assert surgery.measurement == measurement
            assert len(surgery.strip) == d
            assert len(surgery.added) == d + 1
            assert {s.pauli for s in surgery.added} == {new_type}
            assert surgery.removed == []
            assert len(surgery.changed) == d - 1
            for site in surgery.changed:  
                assert site.pauli != new_type
                assert len(old[site]) == 2 and len(surgery.after.checks()[site]) == 4


def test_split_gives_back_the_original_layout():
    for a, b in (("p0", "p1"), ("p0", "p2")):
        layout = grid(3)
        merge = layout.merge(a, b)
        split = merge.after.split(f"{a}+{b}")
        assert split.after.checks() == layout.checks()
        assert split.after.qubits() == layout.qubits()
        assert split.removed == merge.added


def test_only_neighbours_merge():
    with pytest.raises(ValueError):
        grid(3).merge("p0", "p3")


def test_merged_patches_merge_again():
    layout = grid(3).merge("p0", "p1").after
    layout = layout.merge("p2", "p3").after
    surgery = layout.merge("p0+p1", "p2+p3")
    assert surgery.measurement == "XX"
    assert surgery.after.code().k == 1


def test_draw_marks_the_strip_and_new_checks():
    picture = grid(3).merge("p0", "p1").draw()
    assert picture.count("+") == 3
    assert picture.count("z") == 4
