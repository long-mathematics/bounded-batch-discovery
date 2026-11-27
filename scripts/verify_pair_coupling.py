#!/usr/bin/env python3
"""Finite checks for direct pair-clock indicator switching. Requires NumPy and SciPy.
Poisson tail budgets cover omitted terms. All reported values are rounded numerically.
"""
import argparse
import hashlib
import itertools
import json
import math
import platform
from importlib.metadata import version
from pathlib import Path
from collections import defaultdict

import numpy as np
from scipy.stats import norm, poisson


def ensure(condition, label):
    if not condition:
        raise AssertionError(label)


def capped_law(mean, cap):
    return list(poisson.pmf(np.arange(cap), mean)) + [poisson.sf(cap - 1, mean)]


def capped_add(x, y, caps):
    return tuple(min(cap, a + b) for cap, a, b in zip(caps, x, y))


def clock_distribution(edges, caps):
    law = {tuple(0 for _ in caps): 1.0}
    for edge, rate in edges:
        cap = max(caps[j] for j in edge)
        new = defaultdict(float)
        for n, weight in law.items():
            for k, p in enumerate(capped_law(rate, cap)):
                target = tuple(min(c, a + k) if j in edge else a
                               for j, (a, c) in enumerate(zip(n, caps)))
                new[target] += weight * p
        law = new
    return law


def score(n, thresholds, weights):
    return sum(w for x, r, w in zip(n, thresholds, weights) if x >= r)


def marked_laws(edges, i, cap, limit):
    """Neighbor counts after exactly d incident marks for every d up to limit."""
    incident = [(e, a) for e, a in edges if i in e]
    mean = sum(a for _, a in incident)
    increments = []
    for edge, a in incident:
        inc = tuple(int(j in edge and j != i) for j in range(len(cap)))
        increments.append((inc, a / mean))
    laws = [{(0,) * len(cap): 1.0}]
    for _ in range(limit):
        next_law = defaultdict(float)
        for n, p in laws[-1].items():
            for delta, mark_prob in increments:
                next_law[capped_add(n, delta, cap)] += p * mark_prob
        laws.append(next_law)
    return mean, laws


def check_switch(edges, thresholds, weights, limit=19):
    if limit < max(thresholds):
        raise ValueError("Upper Poisson cutoff must reach every threshold.")
    law = clock_distribution(edges, thresholds)
    ylaw = defaultdict(float)
    for n, p in law.items():
        ylaw[score(n, thresholds, weights)] += p
    mean = sum(y * p for y, p in ylaw.items())
    variance = sum((y - mean) ** 2 * p for y, p in ylaw.items())
    v0 = 0.0
    batch_weight = max(sum(weights[i] for i in e) for e, _ in edges)
    cut = sorted(ylaw)[len(ylaw) // 2]
    exact_cov = sum((y - mean) * (y <= cut) * p for y, p in ylaw.items())
    first = second = cov = tail = 0.0
    grid = np.array(sorted(ylaw))
    intervals = defaultdict(float)
    anchored = defaultdict(float)
    for i, (r, weight) in enumerate(zip(thresholds, weights)):
        rate, seq = marked_laws(edges, i, thresholds, limit)
        pmf = poisson.pmf(np.arange(limit + 1), rate)
        p_i = poisson.sf(r - 1, rate)
        v0 += weight ** 2 * p_i * (1 - p_i)
        outside = clock_distribution([(e, a) for e, a in edges if i not in e], thresholds)
        # Unconditioned p_k*p_l on k<r<=l already includes the factor p_i*q_i.
        for k in range(r):
            lower_law = defaultdict(float)
            for n, px in outside.items():
                for counts, py in seq[k].items():
                    lower_law[capped_add(n, counts, thresholds)] += px * py
            # Mix the law of the extra marks in advance; their order remains iid.
            extra = defaultdict(float)
            for l in range(r, limit + 1):
                for counts, p in seq[l - k].items():
                    extra[counts] += pmf[k] * pmf[l] * p
            for n, prob in lower_law.items():
                lo = score(n, thresholds, weights)
                for counts, q in extra.items():
                    hi = weight + score(capped_add(n, counts, thresholds), thresholds, weights)
                    mass = weight * prob * q
                    ensure(hi >= lo, "monotone switch")
                    first += mass * (hi - lo)
                    second += mass * (hi - lo) ** 2
                    length_weight = mass * (hi - lo)
                    intervals[lo, hi] += length_weight
                    anchored[lo] += length_weight * (1 - p_i)
                    anchored[hi] += length_weight * p_i
                    cov += mass * (int(hi <= cut) - int(lo <= cut))
        tail += weight * poisson.cdf(r - 1, rate) * poisson.sf(limit, rate)
    total_weight = sum(weights)
    triangle = any(all(any(set(e) == set(pair) for e, _ in edges)
                       for pair in itertools.combinations(vertices, 2))
                   for vertices in itertools.combinations(range(len(weights)), 3))
    upper = batch_weight * (1.5 * variance + (1.5 if triangle else .5) * v0)
    ensure(abs(first - variance) < 2e-12 + tail * total_weight, "expected switch equals variance")
    ensure(abs(cov - exact_cov) < 2e-12 + tail, "step-test identity")
    ensure(second + tail * total_weight ** 2 < upper + 2e-12, "quadratic-cost theorem")
    crossing = np.array([sum(mass for (lo, hi), mass in intervals.items() if lo <= x < hi)
                         for x in grid])
    exact_probs = np.array([ylaw[y] for y in grid])
    cdf = np.cumsum(exact_probs)
    bias = np.cumsum([anchored[y] for y in grid]) - variance * cdf
    zero_bias_cdf = np.array([-sum((y - mean) * max(x - y, 0) * p for y, p in ylaw.items())
                              / variance for x in grid])
    zero_bias_gap = max(abs(cdf - zero_bias_cdf).max(),
                        abs(cdf - exact_probs - zero_bias_cdf).max())
    normal_cdf = norm.cdf((grid - mean) / math.sqrt(variance))
    dk = max(abs(cdf - normal_cdf).max(), abs(cdf - exact_probs - normal_cdf).max())
    # The omitted unnormalized interval and anchor mass is at most tail*total_weight.
    # Either CDF curve can change by that amount, including at its largest jump.
    omitted_mass = tail * total_weight
    criterion = 2 * (abs(bias).max() + crossing.max() + 2 * omitted_mass) / variance
    uniform_cdf = np.array([
        sum(mass * np.clip((x - lo) / (hi - lo), 0, 1)
            for (lo, hi), mass in intervals.items() if hi > lo) / variance
        for x in grid
    ])
    cdf_error = abs(uniform_cdf - zero_bias_cdf).max()
    ensure(cdf_error <= omitted_mass / variance + 2e-12, "uniform-point zero-bias cdf")
    ensure(abs(np.dot(crossing[:-1], np.diff(grid)) - second) < 2e-12,
           "crossing integral equals quadratic cost")
    ensure(dk <= 2 * zero_bias_gap + 2e-12 and 2 * zero_bias_gap <= criterion + 2e-12,
           "cdf criterion")
    return {"contains_triangle": triangle, "upper_poisson_cutoff": limit,
            "dK": dk, "cdf_criterion_upper": criterion,
            "truncated_crossing_sup": crossing.max(), "truncated_signed_anchor_sup": abs(bias).max(),
            "omitted_length_weight_upper": omitted_mass, "uniform_point_cdf_error": cdf_error,
            "variance": variance, "first_cost": first,
            "step_test_error": abs(cov - exact_cov), "second_cost": second,
            "second_cost_upper_bound": upper,
            "second_cost_truncation_error_upper": tail * total_weight ** 2}


def check_shared_clock():
    # Both threshold indicators switch together. The original law is W*Bernoulli(p)
    # and the zero-bias law is uniform on [0,W], even when the threshold is large.
    r, mean, weights = 24, 24., [.75, 1.25]
    load = sum(weights)
    v = poisson.sf(r - 1, mean) * poisson.cdf(r - 1, mean)
    result = check_switch([((0, 1), mean)], [r, r], weights, limit=85)
    ensure(abs(result["variance"] - load ** 2 * v) < 2e-12, "shared-clock variance")
    ensure(abs(result["second_cost"] - load ** 3 * v) < 2e-12, "shared-clock squared cost")
    ensure(result["truncated_signed_anchor_sup"] < 2e-12, "shared-clock anchor bias")
    return result


def check_cubic():
    edges = [((0, 1), .51), ((1, 2), .38), ((0, 2), .67),
             ((2, 3), .43), ((0,), .24)]
    r = [2, 1, 2, 1]
    h = [1, .5, .75, .25]
    rates = np.array([sum(a for e, a in edges if j in e) for j in range(4)])
    p = poisson.sf(np.array(r) - 1, rates)
    caps = [max(r[j] for j in e) for e, _ in edges]
    marginal = [capped_law(a, cap) for (_, a), cap in zip(edges, caps)]
    configs = np.array(list(itertools.product(*(range(c + 1) for c in caps))))
    prob = np.prod([np.array(marginal[d])[configs[:, d]] for d in range(len(edges))], axis=0)
    touched = [[d for d, (e, _) in enumerate(edges) if j in e] for j in range(4)]
    x = np.array([h[j] * ((configs[:, touched[j]].sum(axis=1) >= r[j]) - p[j])
                  for j in range(4)])

    def projected(j, ids):
        count = configs[:, ids].sum(axis=1)
        remaining = max(0.0, rates[j] - sum(edges[d][1] for d in ids))
        # A clipped count at the cap already guarantees the threshold.
        return h[j] * (poisson.sf(r[j] - count - 1, remaining) - p[j])

    index = {tuple(e): d for d, (e, _) in enumerate(edges) if len(e) == 2}
    first = {}
    second = {}
    for j in range(4):
        adj = [i for i in range(4) if tuple(sorted((i, j))) in index]
        for i in adj:
            first[j, i] = projected(j, [index[tuple(sorted((j, i)))]])
        for i, k in itertools.permutations(adj, 2):
            ids = [index[tuple(sorted((j, i)))], index[tuple(sorted((j, k)))]]
            second[j, i, k] = projected(j, ids) - first[j, i] - first[j, k]

    def expectation(values):
        return float(prob @ values)

    repeat = sum(expectation(x[i] ** 3) for i in range(4))
    repeat += 3 * sum(expectation(x[i] ** 2 * x[j]) for i, j in itertools.permutations(range(4), 2))
    paths = triangles = 0.0
    for i, j, k in itertools.permutations(range(4), 3):
        if (i, j) in first and (k, j) in first:
            paths += 3 * expectation(first[i, j] * x[j] * first[k, j])
        if (i, j, k) in second and (j, k, i) in second and (k, i, j) in second:
            triangles += expectation(second[i, j, k] * second[j, k, i] * second[k, i, j])
    direct = expectation(x.sum(axis=0) ** 3)
    error = abs(direct - repeat - paths - triangles)
    ensure(error < 2e-13 and abs(triangles) > 1e-10, "triangle and path cubic identity")
    return {"cubic": direct, "repeated": repeat, "paths": paths, "triangles": triangles,
            "error": error}


def forest_checks():
    components = []
    for rates, thresholds in [([.8, .37, 1.12], [2, 1, 2, 1]),
                              ([1.2, .86, .51, 1.3], [4, 2, 1, 1, 2]),
                              ([2.1, 1.55], [5, 2, 2])]:
        edges = [((0, j), a) for j, a in enumerate(rates, 1)]
        law = clock_distribution(edges, thresholds)
        dist = np.zeros(len(thresholds) + 1)
        for n, p in law.items():
            value = sum(x >= r for x, r in zip(n, thresholds))
            dist[value] += p
        components.append(dist)
    rows = []
    total = np.ones(1)
    for b in range(1, 61):
        total = np.convolve(total, components[(b - 1) % len(components)])
        if b not in (1, 3, 12, 30, 60):
            continue
        grid = np.arange(len(total))
        mean = total @ grid
        sd = math.sqrt(total @ (grid - mean) ** 2)
        normal = norm.cdf((grid - mean) / sd)
        cdf = np.cumsum(total)
        distance = max(abs(cdf - normal).max(), abs(cdf - total - normal).max())
        rows.append({"components": b, "variance": sd ** 2, "dK": distance,
                     "dK_sqrtV_over_W": distance * sd / 2})
        ensure(abs(total.sum() - 1) < 1e-12, "forest probabilities")
    return rows


def projected_mass_check():
    rows = []
    for r in (16, 64, 256, 1024):
        m = r ** 4
        rate = r / m
        p = -math.expm1(-rate)
        log_clean = poisson.logpmf(r - 1, r) + sum(math.log1p(-k / m) for k in range(r - 1))
        clean = math.exp(log_clean)
        repeated_upper = m * rate ** 2 / 2
        conditional_lower = 2 * (m - r + 1) * p * (1 - p) * clean / (clean + repeated_upper)
        variance_upper = (math.sqrt(m * p * (1 - p)) + .5) ** 2
        ensure(conditional_lower > variance_upper, "clean-atom separation")
        projected_lower = clean * (conditional_lower - variance_upper) ** 2
        rows.append({"hub_threshold": r, "projected_variance_lower_over_r_3_2":
                     projected_lower / r ** 1.5})
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, help="Write the JSON report to this path (default: stdout).")
    args = parser.parse_args()
    report = {
        "switch": [
            check_switch([((0, 1), .55), ((1, 2), .72), ((1, 3), .43), ((2,), .28)],
                         [2, 3, 1, 2], [.5, 1, .75, .25]),
            check_switch([((0, 1), .51), ((1, 2), .38), ((0, 2), .67), ((2, 3), .43),
                          ((0,), .24)], [2, 1, 2, 1], [1, .5, .75, .25])],
        "large_threshold_shared_clock": check_shared_clock(),
        "short_cutoff": check_switch([((0, 1), 1.6), ((0, 2), .9)],
                                     [3, 1, 2], [.5, 1, .25], limit=5),
        "cubic": check_cubic(), "heterogeneous_forest": forest_checks(),
        "projected_switch_mass": projected_mass_check()}
    report["status"] = "passed"
    report["metadata"] = {
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "python_version": platform.python_version(),
        "package_versions": {name: version(name) for name in ('numpy', 'scipy')},
    }
    output = json.dumps(report, indent=2, allow_nan=False)
    if args.output is None:
        print(output)
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(output + "\n", encoding="utf-8")
        print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
