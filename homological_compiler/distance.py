from collections import deque

def distance(layout, pauli="Z"):
    cx = layout.complex() if pauli == "Z" else layout.dual_complex()
    conjugate = "X" if pauli == "Z" else "Z"

    graph = {"R": []}
    for v in cx.cells[0]:
        graph[v] = []
    for edge in cx.cells[1]:
        ends = list(cx.boundary(1, [edge])) + ["R", "R"]  
        u, v = ends[0], ends[1]
        graph[u].append((v, edge))
        graph[v].append((u, edge))

    best = None
    for patch in layout.patches:
        cut = patch.logical(conjugate)
        starts = {u for u in graph for v, edge in graph[u] if edge in cut}
        for start in starts:
            depth = {(start, 0): 0}
            queue = deque([(start, 0)])
            while queue:
                v, layer = queue.popleft()
                for w, edge in graph[v]:
                    step = (w, layer ^ (edge in cut))
                    if step not in depth:
                        depth[step] = depth[v, layer] + 1
                        queue.append(step)
            if (start, 1) in depth and (best is None or depth[start, 1] < best):
                best = depth[start, 1]
    return best
