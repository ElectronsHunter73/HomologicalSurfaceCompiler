import numpy as np
from .css_code import CSSCode
from .geometry import Plaquette, Qubit, rectangle
from .gf2 import GF2Matrix


class Patch:
    def __init__(self,name,row,col,rows,cols,parts=()):
        self.name=name
        self.row,self.col,self.rows,self.cols=row,col,rows,cols
        self.qubits=rectangle(row,col,rows,cols)
        self.parts=parts

    def checks(self):
        checks={}
        for r in range(self.row-1,self.rows+self.row):
            for c in range(self.col-1,self.col+self.cols):
                site=Plaquette(r,c)
                support=frozenset(q for q in site.corners if q in self.qubits)
                if (len(support)==4):
                    checks[site]=support
                elif (len(support)==2):
                    a,b=support
                    vertical=a.c==b.c
                    if(site.pauli=="X")==vertical:
                        checks[site]=support
        return checks
    
    def logical(self,pauli):
        if pauli=="Z":
            return {Qubit(r,self.col) for r  in range(self.row,self.row+self.rows)}
        return {Qubit(self.row,c) for c in range(self.col,self.col+self.cols)}
    

def support_matrix(supports,qubits):
    column={q:i for i,q in enumerate(qubits)}
    m=np.zeros((len(supports),len(qubits)),dtype=np.uint8)
    for i,support in enumerate(supports):
        for q in support:
            m[i,column[q]]=1
    return GF2Matrix(m)

def css_code(checks,qubits):
    hx=support_matrix([checks[s] for s in checks if s.pauli=="X"],qubits)
    hz=support_matrix([checks[s] for s in checks if s.pauli=="Z"],qubits)
    return CSSCode(hx,hz)

def rotated_surface_code(d):
    patch=Patch("patch",0,0,d,d)
    return css_code(patch.checks(),sorted(patch.qubits))
