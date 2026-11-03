import numpy as np
from .css_code import CSSCode
from .gf2 import GF2Matrix

class CellComplex:

    def __init__(self, vertices, edges, faces, d1, d2):
            self.cells = [list(vertices), list(edges), list(faces)]
            self.d = {1: d1, 2: d2}
            self.index = [{cell: i for i, cell in enumerate(cells)} for cells in self.cells]
            if not (d1 @ d2).is_zero():
                raise ValueError("Boundary maps must satisfy d1 @ d2 = 0")

    def chain(self,k,cells):
        vector=np.zeros((len(self.cells[k]),1),dtype=np.uint8)
        for cell in cells:
            vector[self.index[k][cell],0] ^=1
        return GF2Matrix(vector)
    
    def support(self,k,vector):
        return {self.cells[k][i] for i in np.flatnonzero(vector)}

    def cycles(self,k):
        return self.d[k].kernel_basis()

    def boundary(self,k,cells):
        vector= self.d[k] @ self.chain(k,cells)
        return self.support(k-1,vector._data)

    def betti(self,k):
        out_rank=self.d[k].rank() if k in self.d else 0
        in_rank=self.d[k+1].rank() if k+1 in self.d else 0
        return len(self.cells[k]) - out_rank - in_rank

    def dual(self):
        vertices,edges,faces=self.cells
        return CellComplex(faces,edges,vertices,self.d[2].T,self.d[1].T)

    def css_code(self):
        return CSSCode(self.d[1],self.d[2].T)


    def is_subcomplex_of(self,other):
        for k in range(3):
            if any(cell not in other.cells[k] for cell in self.cells[k]):
                return False
        for k in (1,2):
            if any(self.boundary(k,[c]) != other.boundary(k,[c]) for c in self.cells[k]):
                return False
        return True

    def relative_to(self, sub):
        keep = [[i for i,cell in enumerate(self.cells[k]) if cell not in sub.cells[k]] for k in range(3)]
        cells=[[self.cells[k][i] for i in keep[k]] for k in range(3)]
        d1 = GF2Matrix(self.d[1]._data[keep[0], :][:, keep[1]])
        d2 = GF2Matrix(self.d[2]._data[keep[1], :][:, keep[2]])
        return CellComplex(cells[0], cells[1], cells[2], d1, d2)
