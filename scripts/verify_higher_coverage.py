#!/usr/bin/env python3
"""Deterministic finite checks for An optimal variance rate for Poissonized higher coverage.

Run with Python 3, NumPy, SciPy and mpmath. Output is JSON. No random seeds or network.
The threshold distribution is exact up to arithmetic rounding after clock clipping.
Mass checks use a separate larger cutoff with a certified expectation error budget.
"""
import argparse
import hashlib
import itertools
import json
import math
import platform
from importlib.metadata import version
from pathlib import Path
from functools import cache

import mpmath as mp
import numpy as np
from scipy.integrate import quad
from scipy.stats import binom, poisson


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def clipped_poisson(mean, cap):
    """Law of min(X,cap), with the complete Poisson upper tail in the last bin."""
    return np.r_[poisson.pmf(np.arange(cap), mean), poisson.sf(cap - 1, mean)]


def threshold_law(edges, thresholds):
    """Sparse convolution: an exact finite state space for all threshold indicators."""
    states = {tuple(0 for _ in thresholds): 1.0}
    for edge, mean in edges:
        cap = max(thresholds[i] for i in edge)
        new_states = {}
        for state, probability in states.items():
            for count, mass in enumerate(clipped_poisson(mean, cap)):
                target = tuple(min(r, state[i] + (count if i in edge else 0))
                               for i, r in enumerate(thresholds))
                new_states[target] = new_states.get(target, 0.0) + probability * mass
        states = new_states
    return states


def matching_deficit(edges, selected):
    intensity = sum(mean for edge, mean in edges if selected.intersection(edge))
    boundary = {i: sum(mean for edge, mean in edges if selected.intersection(edge) == {i})
                for i in selected}
    pairs = {tuple(sorted(edge)): mean for edge, mean in edges if len(edge) == 2}

    @cache
    def weight(vertices):
        if not vertices:
            return 1.0
        first, *rest = vertices
        result = (1.0 + boundary[first]) * weight(tuple(rest))
        for j in rest:
            rate = pairs.get((first, j), 0.0)
            result += rate * weight(tuple(k for k in rest if k != j))
        return result

    return math.exp(-intensity) * weight(tuple(sorted(selected)))


def check_matching():
    edges = [((0,), .13), ((2,), .21), ((0, 1), .31), ((1, 2), .49),
             ((0, 2), .17), ((2, 3), .23), ((0, 3), .11)]
    states = threshold_law(edges, [2] * 4)
    errors = []
    for selection in itertools.product((False, True), repeat=4):
        ids = {i for i, chosen in enumerate(selection) if chosen}
        exact = sum(p for state, p in states.items() if all(state[i] < 2 for i in ids))
        errors.append(abs(exact - matching_deficit(edges, ids)))
    require(max(errors) < 2e-14, "matching expansion")
    return {"coefficients": len(errors), "max_absolute_error": max(errors)}


def stein_kernel(mean, threshold, counts):
    """From lambda*p_n*H(n) = Cov(1_{X>=r},1_{X>n})."""
    probability = poisson.sf(threshold - 1, mean)
    deficit = poisson.cdf(threshold - 1, mean)
    cov = np.where(counts < threshold, probability * poisson.cdf(counts, mean),
                   deficit * poisson.sf(counts, mean))
    return cov / (mean * poisson.pmf(counts, mean))


def check_mehler():
    errors = []
    for mean, threshold in [(0.41, 2), (2.3, 4), (5.0, 3)]:
        for n in range(8):
            successes = np.arange(min(n, threshold - 1) + 1)

            def integrand(s):
                factors = binom.pmf(successes, n, s)
                return np.dot(factors, poisson.pmf(threshold - 1 - successes, (1 - s) * mean))

            integrated = quad(integrand, 0, 1, epsabs=5e-14)[0]
            errors.append(abs(integrated - stein_kernel(mean, threshold, np.array(n))))
    require(max(errors) < 1e-13, "conditional Mehler inversion")
    return {"cases": len(errors), "max_absolute_error": float(max(errors))}


def check_stein(edges, thresholds, h, cap=17):
    num_edges, species = len(edges), len(h)
    shape = (cap + 1,) * num_edges
    counts = [np.zeros(shape, dtype=int) for _ in h]
    probability = np.ones(shape)
    marginal = np.zeros(species)
    for axis, (edge, mean) in enumerate(edges):
        axis_shape = [1] * num_edges
        axis_shape[axis] = cap + 1
        clock = np.arange(cap + 1).reshape(axis_shape)
        probability *= clipped_poisson(mean, cap).reshape(axis_shape)
        for i in edge:
            counts[i] += clock
            marginal[i] += mean

    def expect(x):
        return np.sum(probability * x)

    y = sum(w * (n >= r) for w, n, r in zip(h, counts, thresholds))
    pivotal = [n == r - 1 for n, r in zip(counts, thresholds)]
    kernel = [stein_kernel(m, r, np.arange(n.max() + 1))[n]
              for m, r, n in zip(marginal, thresholds, counts)]
    variance = expect((y - expect(y)) ** 2)
    mass = np.zeros(shape)
    energy = 0.0
    cut = float(np.median(np.unique(y)))
    test = y <= cut
    cov = expect((y - expect(y)) * test)
    cov_rhs = 0.0
    for edge, rate in edges:
        jump = sum(h[i] * pivotal[i] for i in edge)
        inverse = sum(h[i] * kernel[i] for i in edge)
        mass += rate * jump * inverse
        energy += rate * expect(inverse ** 2)
        # Use exactly representable dyadic weights below for the cdf test.
        cov_rhs += rate * expect(inverse * ((y + jump <= cut).astype(float) - test))
    load = max(sum(h[i] for i in edge) for edge, _ in edges)
    total_rate = sum(rate for _, rate in edges)
    coupling_error = sum(poisson.sf(cap, rate) for _, rate in edges)
    # T and sum a H_C^2 are both in [0, total_rate*load^2].
    maximum_mass = total_rate * load ** 2
    error_budget = coupling_error * maximum_mass
    ell = 1 + sum(poisson.pmf(k, k) for k in range(1, max(thresholds)))
    b = max(thresholds) * poisson.pmf(max(thresholds), max(thresholds))
    mass_var = expect((mass - expect(mass)) ** 2)
    upper = (4 * math.sqrt(ell) + math.sqrt(2 * b)) ** 2 * load ** 2 * variance
    require(abs(expect(mass) - variance) <= error_budget + 3e-13, "mean Stein mass")
    require(energy <= variance + error_budget + 3e-13, "inverse-gradient energy")
    require(abs(cov - cov_rhs) <= 2 * total_rate * load * coupling_error + 3e-13,
            "indicator-test covariance")
    require(mass_var + 3 * maximum_mass ** 2 * coupling_error <= upper + 3e-13,
            "Stein variance bound")
    return {"species": species, "batch_weight": load, "actual_variance": float(variance),
            "mean_mass": float(expect(mass)), "energy": float(energy),
            "mass_variance": float(mass_var), "mass_variance_upper": float(upper),
            "indicator_covariance_error": float(abs(cov - cov_rhs)),
            "omitted_clock_event_probability_upper": float(coupling_error),
            "mass_expectation_truncation_error_upper": float(error_budget)}


def star_fourth(n):
    """Exact finite formulas at high precision, independently of the Gaussian proof.

    E[(S)_j H] obtained by selecting j distinct occupied leaves and
    inclusion-exclusion on which of those j independent clocks are zero.
    """
    def hub_with_zeros(count):
        mean = n - count
        return mp.exp(-count) * mp.gammainc(n, 0, mean, regularized=True)

    def factorial_moment(order, hub):
        falling = math.prod(range(n - order + 1, n + 1))
        if not hub:
            return falling * (1 - mp.exp(-1)) ** order
        return falling * sum((-1) ** i * math.comb(order, i) * hub_with_zeros(i)
                             for i in range(order + 1))

    @cache
    def leaf_moment(order, hub):
        return sum(mp.stirling2(order, j) * factorial_moment(j, hub)
                   for j in range(order + 1))

    raw = [mp.mpf(1)]
    for k in range(1, 5):
        raw.append(leaf_moment(k, False) +
                   sum(math.comb(k, j) * leaf_moment(j, True) for j in range(k)))
    cumulants = [mp.mpf(0)]
    for k in range(1, 5):
        correction = sum(math.comb(k - 1, j - 1) * cumulants[j] * raw[k - j]
                         for j in range(1, k))
        cumulants.append(raw[k] - correction)
    return raw, cumulants


def check_stars():
    mp.mp.dps = 65
    n = 5
    edges = [((0, j), 1.0) for j in range(1, n + 1)]
    states = threshold_law(edges, [n] + [1] * n)
    raw, _ = star_fourth(n)
    for j in range(1, 5):
        independent_check = sum(p * (int(x[0] >= n) + sum(x[1:])) ** j for x, p in states.items())
        require(abs(float(raw[j]) - independent_check) < 1e-11, "star raw moments")
    limit = -4 * mp.exp(-3) / mp.sqrt(2 * mp.pi)
    rows = []
    for n in (32, 128, 512, 2048, 8192):
        _, cumulants = star_fourth(n)
        rows.append({"n": n, "variance_over_n": float(cumulants[2] / n),
                     "fourth_cumulant_over_n_3_2": float(cumulants[4] / mp.mpf(n) ** 1.5)})
    return {"raw_moments_checked_by_clipping": 4, "fourth_cumulant_limit": float(limit),
            "finite_values": rows}


def check_complex_pins():
    """Compute every prefix conditional density for the clocks touching vertex 0."""
    edges = [((0, 1), .48), ((0, 2), .23), ((0, 3), .14),
             ((1, 2), .31), ((2, 3), .52), ((1,), .13)]
    h = np.array([.75, .5, .25, 1.])
    incident = 3
    laws = [clipped_poisson(a, 2) for _, a in edges]
    load = max(sum(h[i] for i in edge) for edge, _ in edges)
    z = (.5 + .8j) / (64 * load * math.sqrt(2))
    s = abs(z) * load
    u = np.expm1(-z * h)
    masses = {}
    deficit = 0j
    for counts in itertools.product(range(3), repeat=len(edges)):
        n = [sum(x for (edge, _), x in zip(edges, counts) if i in edge) for i in range(4)]
        base = math.prod(p[x] for p, x in zip(laws, counts))
        tilted = base * np.prod([1 + u[i] * (n[i] < 2) for i in (1, 2, 3)])
        for length in range(incident + 1):
            prefix = counts[:length]
            masses[prefix] = masses.get(prefix, 0j) + tilted
        if n[0] < 2:
            deficit += tilted
    maximum_ratio = 0.0
    prefixes = 0
    for length in range(incident):
        law = laws[length]
        for prefix in itertools.product(range(3), repeat=length):
            likelihood = np.array([masses[prefix + (x,)] / (masses[prefix] * law[x])
                                   for x in range(3)])
            require(abs(np.dot(law, likelihood) - 1) < 1e-13, "pin normalization")
            excess = np.dot(law, abs(likelihood - 1) ** 2)
            upper = (14 * s) ** 2 * (1 - law[0])
            require(excess <= upper, "pin squared density")
            maximum_ratio = max(maximum_ratio, excess / upper)
            prefixes += 1
    probability = deficit / masses[()]
    require(abs(probability) <= 2, "complex deficit probability")
    return {"prefixes": prefixes, "largest_ratio_to_second_moment_allowance": float(maximum_ratio),
            "tilted_deficit_real": probability.real, "tilted_deficit_imag": probability.imag}


def check_projected_mass():
    # Lower bound from the single clean atom, with an upper bound on the variance
    # of Y. O(r) scalar operations even though the star has r^4 leaves.
    mp.mp.dps = 65
    rows = []
    for r in (16, 64, 256, 1024):
        m = r ** 4
        a = mp.mpf(r) / m
        occupied = -mp.expm1(-a)
        clean = mp.exp(-r) * mp.mpf(r) ** (r - 1) / mp.factorial(r - 1)
        clean *= mp.fprod(1 - mp.mpf(i) / m for i in range(r - 1))
        collision_upper = m * a ** 2 / 2
        conditional_lower = 2 * (m - r + 1) * occupied * clean / (clean + collision_upper)
        variance_upper = (mp.sqrt(m * occupied * (1 - occupied)) + mp.mpf('.5')) ** 2
        require(conditional_lower > variance_upper, "separation of single projected atom")
        lower = clean * (conditional_lower - variance_upper) ** 2
        rows.append({"r": r, "projected_variance_lower_over_r_3_2": float(lower / mp.mpf(r) ** 1.5)})
    return {"asymptotic_lower_constant": float(1 / mp.sqrt(2 * mp.pi)), "finite_bounds": rows}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, help="Write the JSON report to this path (default: stdout).")
    args = parser.parse_args()
    report = {"threshold_two_matching": check_matching(), "Mehler": check_mehler(),
              "Stein": [
                  check_stein([((0,), .45), ((0, 1), .11), ((1,), .23)],
                              [2, 2], [.75, .5]),
                  check_stein([((0, 1, 2), .06), ((0, 1), .05), ((0, 2), .04),
                               ((0,), .85)], [3, 2, 1], [1, .5, .25])],
              "complex_pins": check_complex_pins(), "star": check_stars(),
              "projected_mass": check_projected_mass()}
    report["status"] = "passed"
    report["metadata"] = {
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "python_version": platform.python_version(),
        "package_versions": {name: version(name) for name in ('numpy', 'scipy', 'mpmath')},
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
