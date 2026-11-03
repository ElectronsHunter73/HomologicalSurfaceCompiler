from .chain_complex import CellComplex
from .geometry import bounding_box, rectangle
from .rotated_surface import Patch, css_code, support_matrix

def strip_between(a,b):
    a,b = sorted([a,b], key=lambda p:(p.row,p.col))
    if (a.row,a.rows) == (b.row,b.rows) and b.col == a.col + a.cols +1:
        return rectangle(a.row,a.col+a.cols,a.rows,1),"ZZ"
    if (a.col,a.cols) == (b.col,b.cols) and b.row == a.row + a.rows +1:
        return rectangle(a.row+a.rows,a.col,1,a.cols), "XX"
    return None, None 

class SurfaceLayout:

    def __init__(self,patches=()):
        self.patches=list(patches)

    def patch(self,name):
        for p in self.patches:
            if p.name==name:
                return p

    def add_patch(self,name,d,next_to=None,side=None):
        row,col=0,0
        if next_to is not None:
            other=self.patch(next_to)
            row,col={
                "left":(other.row,other.col-d-1),
                "right":(other.row,other.col+other.cols+1),
                "top":(other.row-d-1,other.col),
                "bottom":(other.row+other.rows+1,other.col)}[side]
        self.patches.append(Patch(name,row,col,d,d))
    
    def qubits(self):
        return sorted(q for p in self.patches for q in p.qubits)

    def checks(self):
        checks ={}
        for patch in self.patches:
            for site,support in patch.checks().items():
                if site in checks:
                    raise ValueError(f"Two patches use the check site {site}")
                checks[site]=support
        return dict(sorted(checks.items()))

    def vector(self,qubits):
        return support_matrix([qubits],self.qubits()).T

    def code(self):
        return css_code(self.checks(),self.qubits())

    def complex(self):
        checks,code = self.checks(), self.code()
        x_sites = [s for s in checks if s.pauli =="X"]
        z_sites = [s for s in checks if s.pauli =="Z"]
        return CellComplex(x_sites,self.qubits(),z_sites,code.hx,code.hz.T)

    def dual_complex(self):
        return self.complex().dual()

    def neighbours(self):
        pairs=[]
        for i,a in enumerate(self.patches):
            for b in self.patches[i+1:]:
                strip,measurement=strip_between(a,b)
                if strip:
                    pairs.append((a.name,b.name,measurement))
        return pairs

    def merge(self,a,b):
        pa,pb=self.patch(a),self.patch(b)
        strip,measurement=strip_between(pa,pb)
        if not strip:
            raise ValueError(f"{a} and {b} do not face each other across a free strip")
        row,col = min(pa.row, pb.row), min(pa.col, pb.col)
        rows= max(pa.row+pa.rows,pb.row+pb.rows) - row
        cols= max(pa.col+pa.cols,pb.col+pb.cols) - col
        merged=Patch(f"{a}+{b}",row,col,rows,cols,(pa,pb))
        after = SurfaceLayout([merged if p.name == a else p for p in self.patches if p.name != b])
        return Surgery("merge", measurement, self, after, [a, b], [merged.name], strip)

    def split(self,name):
        pa,pb=self.patch(name).parts
        strip,measurement=strip_between(pa,pb)
        after =SurfaceLayout([p for p in self.patches if p.name != name] +[pa,pb])
        return Surgery("split", measurement, self, after, [name], [pa.name, pb.name], strip)

    def draw(self, strip=(), new_checks=()):
        """o = qubit, X/Z = check; strip qubits as +, new checks in lowercase."""
        grid = {}
        for q in self.qubits():
            grid[2 * q.r, 2 * q.c] = "+" if q in strip else "o"
        for s in self.checks():
            grid[2 * s.r + 1, 2 * s.c + 1] = s.pauli.lower() if s in new_checks else s.pauli
        r0, c0, r1, c1 = bounding_box(self.qubits())
        rows = range(2 * r0 - 1, 2 * r1 + 2)
        cols = range(2 * c0 - 1, 2 * c1 + 2)
        return "\n".join(" ".join(grid.get((y, x), " ") for x in cols).rstrip() for y in rows)


class Surgery:
    def __init__(self,kind,measurement,before,after,inputs,outputs,strip):
        self.kind,self.measurement=kind,measurement
        self.before,self.after=before,after
        self.inputs,self.outputs=inputs,outputs
        self.strip=strip
        self.old_checks,self.new_checks=before.checks(),after.checks()
        old,new=self.old_checks,self.new_checks
        self.added = [s for s in new if s not in old]
        self.removed = [s for s in old if s not in new]
        self.changed = [ s for s in new if s in old and old[s] != new[s]]

    def strip_action(self):
        if self.kind == "merge":
            return "prepare |+>" if self.measurement == "ZZ" else "prepare |0>"
        return "measure X" if self.measurement == "ZZ" else "measure Z"

    def summary(self):
        old, new = self.old_checks, self.new_checks
        changed = [f"{s} w{len(old[s])}->w{len(new[s])}" for s in self.changed]
        return "\n".join([
            f"{self.kind} {' + '.join(self.inputs)} -> {', '.join(self.outputs)} ({self.measurement})",
            f"  strip: {len(self.strip)} qubits, {self.strip_action()}",
            f"  added checks:   {self.added}",
            f"  removed checks: {self.removed}",
            f"  changed checks: {changed}",])

    def draw(self):
            if self.kind == "merge":
                return self.after.draw(self.strip, self.added)
            return self.before.draw(self.strip, self.removed)


      
        




    