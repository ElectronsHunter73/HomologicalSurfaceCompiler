from homological_compiler.distance import distance
from homological_compiler.homology import Prediction, gauge_fix
from homological_compiler.layout import SurfaceLayout
from homological_compiler.verification import verify

layout = SurfaceLayout()
layout.add_patch("p0", 3)
layout.add_patch("p1", 3, next_to="p0", side="right")
layout.add_patch("p2", 3, next_to="p0", side="bottom")
layout.add_patch("p3", 3, next_to="p1", side="bottom")

print(layout.draw())
print("\nneighbours:", layout.neighbours())

merge = layout.merge("p0", "p1")
print("\n" + merge.summary())
print(merge.draw())

prediction = Prediction(merge)
print("\n" + prediction.report())
print("GF(2) verification:", verify(prediction) or "all predictions confirmed")

outcomes = {site: 0 for site in merge.added}
outcomes[merge.added[0]] = outcomes[merge.added[1]] = 1 
m, correction = gauge_fix(prediction, outcomes)
print(f"\noutcomes {outcomes}\n  logical result m = {m}, correction X on {sorted(correction)}")

print("\ndistance after the merge: Z", distance(merge.after, "Z"), " X", distance(merge.after, "X"))

split = merge.after.split("p0+p1")
print("\n" + Prediction(split).report())
print("split restores the layout:", split.after.checks() == layout.checks())
