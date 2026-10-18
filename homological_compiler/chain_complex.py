import numpy as np
from .gf2 import GF2Matrix
from .rotated_surface import rotated_surface_code


class CellComplex:
   
    def __init__(self, n_vertices, edges, faces, boundary_vertices=()):
        self.vertices = tuple(range(n_vertices))
        self.edges = tuple(tuple(edge) for edge in edges)
        self.faces = tuple(tuple(face) for face in faces)
        self.boundary_vertices = tuple(boundary_vertices)
        if any(v not in self.vertices for v in self.boundary_vertices):
            raise ValueError("Invalid boundary vertex")
        d1 = np.zeros((len(self.vertices), len(self.edges)),dtype=np.uint8)
        for j, (u, v) in enumerate(self.edges):
            if (u not in self.vertices or v not in self.vertices or u == v):
                raise ValueError("Invalid edge endpoints")
            d1[u, j] = 1
            d1[v, j] = 1
        d2 = np.zeros((len(self.edges), len(self.faces)),dtype=np.uint8)
        for j, support in enumerate(self.faces):
            for edge in support:
                if not (0 <= edge < len(self.edges)):
                    raise ValueError("Invalid face edge")
                d2[edge, j] ^= 1
        self.boundary_1 = GF2Matrix(d1)
        self.boundary_2 = GF2Matrix(d2)
        if not (self.boundary_1 @ self.boundary_2).is_zero():
            raise ValueError("Boundary maps must satisfy d1 @ d2 = 0")
        
    def relative_boundary_maps(self):
        vertices = [v for v in self.vertices if v not in self.boundary_vertices]
        d1 = GF2Matrix(self.boundary_1._data[vertices, :])
        d2 = self.boundary_2
        return d1, d2

    def relative_h1_dimension(self):
        d1, d2 = self.relative_boundary_maps()
        return d1.n_cols - d1.rank() - d2.rank()



def rotated_surface_cell_complex(d):
    code = rotated_surface_code(d)
    top = code.hx.n_rows
    bottom = top + 1
    edges = []
    for qubit in range(code.n):
        checks = np.flatnonzero(code.hx._data[:, qubit])
        if len(checks) == 2:
            edges.append((int(checks[0]), int(checks[1])))
        elif len(checks) == 1:
            row = qubit // d
            if row == 0:
                boundary = top
            elif row == d - 1:
                boundary = bottom
            else:
                raise ValueError("Unexpected boundary qubit")
            edges.append((int(checks[0]), boundary))
        else:
            raise ValueError("Invalid qubit incidence")
    faces = []
    for row in code.hz._data:
        support = np.flatnonzero(row).tolist()
        faces.append(support)
    return CellComplex(n_vertices=code.hx.n_rows + 2,edges=edges,faces=faces,boundary_vertices=(top, bottom))
