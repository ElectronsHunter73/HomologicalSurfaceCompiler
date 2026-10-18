from .chain_complex import CellComplex, rotated_surface_cell_complex
from .css_code import CSSCode

class PatchCollection:
    def __init__(self):
         self.patches = []

    def add_patch(self, d):
        patch = rotated_surface_cell_complex(d)
        self.patches.append(patch)

    def cell_complex(self):
        if not self.patches:
            raise ValueError("No patches to combine")
        edges = []
        faces = []
        boundary_vertices = []
        vertex_offset = 0
        edge_offset = 0

        for patch in self.patches:
            for u, v in patch.edges:
                edges.append((u + vertex_offset,v + vertex_offset))

            for support in patch.faces:
                faces.append([edge + edge_offset for edge in support])

            for vertex in patch.boundary_vertices:
                boundary_vertices.append(vertex + vertex_offset)

            vertex_offset += len(patch.vertices)
            edge_offset += len(patch.edges)

        return CellComplex(n_vertices=vertex_offset,edges=edges,faces=faces,boundary_vertices=boundary_vertices)

    def css_code(self):
        cells = self.cell_complex()
        d1, d2 = cells.relative_boundary_maps()
        return CSSCode(d1, d2.T)
