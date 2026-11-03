import pytest

from homological_compiler.chain_complex import CellComplex
from homological_compiler.gf2 import GF2Matrix
from homological_compiler.layout import SurfaceLayout


def two_patches(d=3):
    layout = SurfaceLayout()
    layout.add_patch("a", d)
    layout.add_patch("b", d, next_to="a", side="right")
    return layout


def test_boundary_of_boundary_must_vanish():
    with pytest.raises(ValueError):
        CellComplex(["v"], ["e"], ["f"], GF2Matrix([[1]]), GF2Matrix([[1]]))


def test_the_complex_is_the_code():
    for d in (3, 5):
        layout = two_patches(d)
        cx, code = layout.complex(), layout.code()
        assert cx.css_code().hx == code.hx and cx.css_code().hz == code.hz
        assert cx.betti(1) == cx.dual().betti(1) == code.k == 2 
        assert cx.betti(0) == 0  


def test_dual_twice_is_the_identity():
    cx = two_patches().complex()
    back = cx.dual().dual()
    assert back.cells == cx.cells and back.d[1] == cx.d[1] and back.d[2] == cx.d[2]


def test_chains_are_sets_of_named_cells():
    layout = two_patches()
    cx = layout.complex()
    assert cx.boundary(1, layout.patch("a").logical("Z")) == set() 
    for site, support in layout.checks().items():
        if site.pauli == "Z":  
            assert cx.boundary(2, [site]) == support


def test_zz_merge_is_an_inclusion_of_primal_complexes():
    surgery = two_patches().merge("a", "b")
    Y, X = surgery.before.complex(), surgery.after.complex()
    assert Y.is_subcomplex_of(X)
    assert not Y.dual().is_subcomplex_of(X.dual())
    Q = X.relative_to(Y)
    assert set(Q.cells[1]) == set(surgery.strip)
    assert set(Q.cells[2]) == set(surgery.added)
    assert Q.cells[0] == []  
    assert Q.betti(2) == 1 and Q.betti(1) == 0  
