from dataclasses import dataclass

@dataclass(frozen=True,order=True)
class Qubit:
    r:int
    c:int

    def __repr__(self):
            return f"q({self.r},{self.c})"

@dataclass(frozen=True,order=True)
class Plaquette:
      r:int
      c:int
      
      @property
      def pauli(self):
            return "X" if (self.r+self.c)%2 ==0 else "Z"
      @property
      def corners(self):
            r,c=self.r,self.c
            return (Qubit(r,c),Qubit(r,c+1),Qubit(r+1,c),Qubit(r+1,c+1))
      
      def __repr__(self) :
              return f"{self.pauli}({self.r},{self.c})"
      
def rectangle(row,col,rows,cols):
      return frozenset(Qubit(r,c) for r in range(row,row+rows) for c in range(col,col+cols))
      
            
def bounding_box(qubits) -> tuple[int, int, int, int]:
    rows = [q.r for q in qubits]
    cols = [q.c for q in qubits]
    return min(rows), min(cols), max(rows), max(cols)