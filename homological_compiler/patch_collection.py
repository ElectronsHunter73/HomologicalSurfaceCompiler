from .chain_complex import CellComplex, rotated_surface_cell_complex
from .css_code import CSSCode
from .rotated_surface import qubit_index

class PatchCollection:
    def __init__(self):
         self.patches = []
         self.distances = []
         self.positions =[]

    def add_patch(self, d, next_to=None, side=None):
        if not self.patches:
            if next_to is not None or side is not None:
                raise ValueError("First patch needs no reference")
            position=(0,0)

        elif next_to is None and side is None:
            position = None

        else:
            if next_to is None or side is None:
                raise ValueError("Specify both next_to and side")

            if not 0 <= next_to <len(self.patches):
                raise ValueError("Invalid reference patch")

            if self.positions[next_to] is None:
                raise ValueError("The reference patch is not placed")

            if d != self.distances[next_to]:
                raise ValueError("Adjancent patches must have equal distance")

            row,col = self.positions[next_to]
            if side == "left":
                position = (row, col - d)

            elif side == "right":
                position = (row, col + d)

            elif side == "top":
                position = (row - d, col)

            elif side == "bottom":
                position = (row + d, col)

            else:
                raise ValueError("Invalid side")

            new_row, new_col= position
            for i, other in enumerate(self.positions):
                if other is None:
                    continue
                other_row, other_col = other
                other_d = self.distances[i]
                overlap = (
                    new_row < other_row + other_d
                    and other_row < new_row + d
                    and new_col < other_col + other_d
                    and other_col < new_col + d
                )
                if overlap:
                    raise ValueError("Patches cannot overlap")
        patch=rotated_surface_cell_complex(d)
        self.patches.append(patch)
        self.distances.append(d)
        self.positions.append(position)
        return len(self.patches)-1

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
    
    def boundary_qubits(self, patch_index, side):
        if not 0 <= patch_index < len(self.patches):
            raise IndexError("Invalid patch index")
        d = self.distances[patch_index]

        if side == "top":
            coords = [(0, col) for col in range(d)]

        elif side == "bottom":
            coords = [(d - 1, col) for col in range(d)]

        elif side == "left":
            coords = [(row, 0) for row in range(d)]

        elif side == "right":
            coords = [(row, d - 1) for row in range(d)]
        else:
            raise ValueError("Invalid boundary side")

        offset = sum(len(patch.edges) for patch in self.patches[:patch_index])

        return tuple(offset + qubit_index(row, col, d) for row, col in coords)
    def facing_sides(self, patch_a, patch_b):

        if patch_a == patch_b:
            raise ValueError("Choose two different patches")

        for i in (patch_a, patch_b):
            if not 0 <= i < len(self.patches):
                raise IndexError("Invalid patch index")

        pos_a = self.positions[patch_a]
        pos_b = self.positions[patch_b]

        if pos_a is None or pos_b is None:
            return None

        d_a = self.distances[patch_a]
        d_b = self.distances[patch_b]

        if d_a != d_b:
            return None

        row_a, col_a = pos_a
        row_b, col_b = pos_b

        if row_a == row_b:
            if col_a + d_a == col_b:
                return ("right", "left")

            if col_b + d_b == col_a:
                return ("left", "right")

        if col_a == col_b:
            if row_a + d_a == row_b:
                return ("bottom", "top")

            if row_b + d_b == row_a:
                return ("top", "bottom")

        return None

    def adjacent_pairs(self):
        pairs = []
        for i in range(len(self.patches)):
            for j in range(i + 1, len(self.patches)):

                sides = self.facing_sides(i, j)

                if sides is not None:
                    pairs.append((i, j, sides[0], sides[1]))

        return pairs
    
    def candidate_interface(self, patch_a, patch_b):
        sides = self.facing_sides(patch_a, patch_b)
        if sides is None:
            return None
        side_a, side_b = sides
        if side_a in ("left", "right"):
            measurement = "ZZ"
        else:
            measurement = "XX"
        qubits_a = self.boundary_qubits(patch_a, side_a)
        qubits_b = self.boundary_qubits(patch_b, side_b)
        pairs = tuple(zip(qubits_a, qubits_b))
        return measurement, pairs

    def touching_stabilizers(self, patch_a, patch_b):
        sides = self.facing_sides(patch_a, patch_b)
        if sides is None:
            return None
        side_a, side_b = sides
        boundary = set(self.boundary_qubits(patch_a, side_a) + self.boundary_qubits(patch_b, side_b))
        code = self.css_code()
        result = {}
        for kind, matrix in (("X", code.hx), ("Z", code.hz)):
            result[kind] = tuple(i for i, row in enumerate(matrix._data) if any(row[q] for q in boundary))
        return result

    def interfaces_commute(self, first, second):
        measurement_a = self.candidate_interface(*first)
        measurement_b = self.candidate_interface(*second)
        if measurement_a is None or measurement_b is None:
            raise ValueError("Both interfaces must be adjacent")
        type_a = measurement_a[0]
        type_b = measurement_b[0]
        if type_a == type_b:
            return True
        shared = len(set(first) & set(second))
        return shared % 2 == 0

