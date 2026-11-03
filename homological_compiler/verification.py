
def verify(prediction):
    s, pauli = prediction.surgery, prediction.pauli
    small, big = prediction.small.code(), prediction.big.code()
    before, after = (small, big) if s.kind == "merge" else (big, small)
    problems = []
    if before.k != prediction.k_before:
        problems.append(f"k before is {before.k}, predicted {prediction.k_before}")
    if after.k != prediction.k_after:
        problems.append(f"k after is {after.k}, predicted {prediction.k_after}")

    checks = prediction.big.checks()
    for label, op, S in zip(prediction.labels, prediction.measured, prediction.generators):
        if pauli == "Z":
            stabilizer = big.is_z_stabilizer(prediction.big.vector(op))
            logical = small.is_z_logical(prediction.small.vector(op))
        else:
            stabilizer = big.is_x_stabilizer(prediction.big.vector(op))
            logical = small.is_x_logical(prediction.small.vector(op))
        if not stabilizer:
            problems.append(f"{label} is not a stabilizer of the larger code")
        if not logical:
            problems.append(f"{label} is not a logical operator of the smaller code")
        product = set()
        for site in S:
            product ^= checks[site]
        if product != op:
            problems.append(f"the checks recorded for {label} do not multiply to it")

    if s.measurement != pauli * 2:
        problems.append(f"geometry says {s.measurement}, homology says {pauli * 2}")
    return problems
