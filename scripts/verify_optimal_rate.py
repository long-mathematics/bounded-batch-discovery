#!/usr/bin/env python3
"""Deterministic finite-model checks for the optimal Poissonized rate.

Requires mpmath. Exhaustive batch-presence enumeration, independent-set sums,
and exact finite-sum formulas are evaluated at 80 decimal digits. Fixed complex
activities and induced subgraphs supplement the proofs; they do not establish
zero-freeness or the Berry--Esseen theorem. No Monte Carlo sampling is used.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
from functools import lru_cache
import hashlib
from importlib.metadata import version
import json
from math import comb, factorial
from pathlib import Path
import platform

import mpmath as mp


if not __debug__:
    raise SystemExit('Verification requires assertions; run Python without -O or -OO.')

PRECISION = 80
TOLERANCE = '1e-60'
TIMES = ('0.05', '1', '8')
MAX_CUMULANT_ORDER = 8


def scalar(value):
    value = Fraction(value)
    return mp.mpf(value.numerator) / value.denominator


def members(mask):
    while mask:
        bit = mask & -mask
        yield bit.bit_length() - 1
        mask ^= bit


def models():
    dense = [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3),
             (0, 1), (1, 2), (2, 3), (0, 3)]
    star = [(0, i) for i in range(1, 5)]
    blocks = [(0, 1), (2, 3, 4), (5,)]
    # Integer rate weights are normalized to batch probabilities below.
    return [
        ('singletons', [(0,), (1,), (2,), (3,)], [1, 3, 6, 10], ['1'] * 4),
        ('disjoint_blocks', blocks, [2, 5, 3], ['1'] * 6),
        ('overlapping_star', star, [7, 13, 30, 50], ['1'] * 5),
        ('dense_triples', dense, list(range(1, 9)), ['1'] * 4),
        ('empty_and_zero_rate', [(), (0,), (1, 2), (0, 1, 3),
                                 (2, 3, 4), (4,), (1, 4)],
         [1, 2, 4, 3, 6, 4, 0], ['1'] * 5),
        ('weighted_star', star, [7, 13, 30, 50],
         ['3/2', '1/4', '2', '0', '3/4']),
        ('weighted_dense', dense, list(range(1, 9)),
         ['5/2', '1/8', '3/8', '5/4']),
        ('weighted_disjoint', blocks, [2, 5, 3],
         ['1/4', '3/4', '3/2', '1/2', '1/8', '7/2']),
        ('weighted_many_small', [tuple(range(6)), (0, 1), (1, 2),
                                 (2, 3), (3, 4), (4, 5), (0, 5)],
         [2, 1, 2, 3, 4, 5, 3], ['1/16', '3/16', '1/8', '0', '1/4', '3/8']),
    ]


class IncidenceGraph:
    """Independent sets are enumerated directly, without ratio recursion."""

    def __init__(self, batches, species_count):
        self.d = species_count
        self.n = species_count + len(batches)
        self.full = (1 << self.n) - 1
        self.neighbors = [0] * self.n
        for j, batch in enumerate(batches, species_count):
            for i in batch:
                self.neighbors[i] |= 1 << j
                self.neighbors[j] |= 1 << i
        self.independent = []
        for mask in range(1 << self.n):
            vertices = tuple(members(mask))
            if all(not (self.neighbors[v] & mask) for v in vertices):
                self.independent.append((mask, vertices))

    def partition_function(self, activities):
        terms = [(mask, mp.fprod(activities[v] for v in vertices))
                 for mask, vertices in self.independent]

        @lru_cache(maxsize=None)
        def partition(mask):
            return mp.fsum(value for subset, value in terms if subset & mask == subset)

        return partition

    def selected_subgraphs(self):
        # Fixed proper induced subgraphs, including isolated-vertex and mixed cases.
        masks = {0, self.full, (1 << self.d) - 1,
                 self.full ^ ((1 << self.d) - 1)}
        masks.update(self.full ^ (1 << v) for v in range(self.n))
        masks.update(sum(1 << v for v in range(self.n) if v % 3 != j)
                     for j in range(3))
        masks.update(sum(1 << v for v in range(self.n) if v % 2 == j)
                     for j in range(2))
        return sorted(masks)


def presence_distribution(batch_masks, intensities):
    """Independent Poisson clocks reduce exactly to Bernoulli presence states."""
    absent = [mp.exp(-a) for a in intensities]
    present = [-mp.expm1(-a) for a in intensities]
    law = {}
    for state in range(1 << len(batch_masks)):
        probability = mp.mpf(1)
        seen = 0
        for j, batch in enumerate(batch_masks):
            if state & (1 << j):
                probability *= present[j]
                seen |= batch
            else:
                probability *= absent[j]
        assert probability >= 0
        law[seen] = law.get(seen, mp.mpf(0)) + probability
    return law


def absence_polynomial(activities, batch_masks, intensities):
    # Inclusion probabilities of arbitrary absent species sets, independently
    # derived from their incident Poisson clocks.
    return mp.fsum(
        mp.fprod(activities[i] for i in members(subset))
        * mp.exp(-mp.fsum(a for batch, a in zip(batch_masks, intensities)
                         if batch & subset))
        for subset in range(1 << len(activities)))


def relative_error(left, right):
    return abs(left - right) / max(mp.mpf(1), abs(left), abs(right))


def check_ratios(graph, partition, activities, radii, totals, maxima):
    tolerance = mp.mpf(TOLERANCE)
    for mask in graph.selected_subgraphs():
        value = partition(mask)
        assert abs(value) > tolerance
        totals['induced_subgraph_checks'] += 1
        totals['proper_induced_subgraph_checks'] += int(mask != graph.full)
        for v in members(mask):
            without_v = mask ^ (1 << v)
            ratio = activities[v] * partition(without_v & ~graph.neighbors[v]) / partition(without_v)
            # Compare the directly enumerated partition ratio to the telescoping
            # product on smaller neighbor-deleted induced subgraphs.
            product = activities[v]
            current = without_v
            for neighbor in members(graph.neighbors[v] & mask):
                after = current ^ (1 << neighbor)
                child = activities[neighbor] * partition(after & ~graph.neighbors[neighbor]) / partition(after)
                assert abs(1 + child) > tolerance
                product /= 1 + child
                current = after
            error = relative_error(ratio, product)
            assert error < tolerance
            maxima['ratio_recursion_relative_error'] = max(maxima['ratio_recursion_relative_error'], error)
            assert relative_error(value, partition(without_v) * (1 + ratio)) < tolerance
            totals['ratio_identity_checks'] += 1
            if v < graph.d:
                assert abs(ratio) <= radii[v] + tolerance
                totals['species_radius_checks'] += 1
                if radii[v] > 0:
                    maxima['species_ratio_over_radius'] = max(
                        maxima['species_ratio_over_radius'], abs(ratio) / radii[v])
            else:
                assert mp.re(ratio) >= -tolerance
                totals['batch_nonnegative_real_checks'] += 1
                maxima['batch_negative_real_violation'] = max(
                    maxima['batch_negative_real_violation'], -mp.re(ratio))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path,
                        default=Path(__file__).resolve().parents[1] / 'reports' / 'optimal_rate_verification.json',
                        help='JSON report path (default: repository reports/optimal_rate_verification.json)')
    args = parser.parse_args()
    mp.mp.dps = PRECISION
    tolerance = mp.mpf(TOLERANCE)
    totals = {key: 0 for key in (
        'finite_models', 'weighted_models', 'distribution_checks', 'mean_identity_checks',
        'variance_identity_checks', 'cumulant_bound_checks', 'complex_partition_identity_checks',
        'absence_expectation_identity_checks', 'weighted_laplace_identity_checks',
        'batch_sector_checks', 'induced_subgraph_checks', 'proper_induced_subgraph_checks',
        'ratio_identity_checks', 'species_radius_checks', 'batch_nonnegative_real_checks')}
    maxima = {key: mp.mpf(0) for key in (
        'normalization_error', 'mean_relative_error', 'variance_relative_error',
        'partition_identity_relative_error', 'absence_expectation_relative_error',
        'weighted_laplace_relative_error', 'ratio_recursion_relative_error',
        'species_ratio_over_radius', 'batch_negative_real_violation', 'cumulant_bound_ratio')}
    specifications, records = [], []
    for name, batches, rate_weights, weight_strings in models():
        d = len(weight_strings)
        assert len(batches) == len(rate_weights) and sum(rate_weights) > 0
        assert all(rate >= 0 for rate in rate_weights)
        assert all(len(set(batch)) == len(batch) and all(0 <= i < d for i in batch)
                   for batch in batches)
        batch_masks = [sum(1 << i for i in batch) for batch in batches]
        species_weights = [scalar(w) for w in weight_strings]
        assert all(w >= 0 for w in species_weights)
        B = max(map(len, batches))
        W = max(mp.fsum(species_weights[i] for i in batch) for batch in batches)
        assert W > 0
        graph = IncidenceGraph(batches, d)
        weighted = any(Fraction(w) != 1 for w in weight_strings)
        totals['finite_models'] += 1
        totals['weighted_models'] += int(weighted)
        specifications.append({'name': name, 'batches': batches, 'rate_weights': rate_weights,
                               'species_weights': weight_strings, 'B': B, 'W': float(W),
                               'independent_sets_enumerated': len(graph.independent)})
        for time_string in TIMES:
            t = scalar(time_string)
            intensities = [t * rate / sum(rate_weights) for rate in rate_weights]
            seen_law = presence_distribution(batch_masks, intensities)
            normalization = abs(mp.fsum(seen_law.values()) - 1)
            assert normalization < tolerance
            maxima['normalization_error'] = max(maxima['normalization_error'], normalization)
            totals['distribution_checks'] += 1
            # Rational scores group states into the exact finite support before
            # central moments and cumulants are evaluated at high precision.
            pmf = {}
            for seen, probability in seen_law.items():
                score = sum((Fraction(weight_strings[i]) for i in members(seen)), Fraction(0))
                pmf[score] = pmf.get(score, mp.mpf(0)) + probability
            mean = mp.fsum(scalar(x) * p for x, p in pmf.items())
            moments = [mp.fsum((scalar(x) - mean)**k * p for x, p in pmf.items())
                       for k in range(MAX_CUMULANT_ORDER + 1)]
            V = moments[2]
            assert V > 0
            absent_i = [mp.exp(-mp.fsum(a for batch, a in zip(batch_masks, intensities)
                                       if batch & (1 << i))) for i in range(d)]
            closed_mean = mp.fsum(w * (1 - q) for w, q in zip(species_weights, absent_i))
            closed_variance = mp.fsum(
                species_weights[i] * species_weights[j]
                * (mp.exp(-mp.fsum(a for batch, a in zip(batch_masks, intensities)
                                   if batch & ((1 << i) | (1 << j)))) - absent_i[i] * absent_i[j])
                for i in range(d) for j in range(d))
            for label, left, right in (('mean', mean, closed_mean), ('variance', V, closed_variance)):
                error = relative_error(left, right)
                assert error < tolerance
                maxima[label + '_relative_error'] = max(maxima[label + '_relative_error'], error)
                totals[label + '_identity_checks'] += 1
            cumulants = [mp.mpf(0)]
            for k in range(1, MAX_CUMULANT_ORDER + 1):
                cumulants.append(moments[k] - mp.fsum(
                    comb(k - 1, j - 1) * cumulants[j] * moments[k - j] for j in range(1, k)))
            assert abs(cumulants[1]) < tolerance and abs(cumulants[2] - V) < tolerance
            ratios = []
            for k in range(3, MAX_CUMULANT_ORDER + 1):
                bound = 2 * factorial(k) * V * (4 * W)**(k - 2)
                ratio = abs(cumulants[k]) / bound
                assert ratio <= 1 + tolerance, (name, time_string, k, ratio)
                ratios.append(float(ratio))
                maxima['cumulant_bound_ratio'] = max(maxima['cumulant_bound_ratio'], ratio)
                totals['cumulant_bound_checks'] += 1
            records.append({'model': name, 'time': float(t), 'pmf_support_size': len(pmf),
                            'mean': float(mean), 'variance': float(V),
                            'cumulant_bound_ratios_orders_3_to_8': ratios})

            radius = mp.mpf('0.95') * mp.sin(mp.pi / (2 * B))
            profiles = [
                ('polydisk_phases', [radius * mp.exp(mp.j * (i + 1) * mp.pi / 7) for i in range(d)],
                 [radius] * d, None),
                ('polydisk_negative', [-radius] * d, [radius] * d, None),
            ]
            for label, scaled_theta in (
                    ('laplace_positive', mp.mpf('0.49')),
                    ('laplace_negative', mp.mpf('-0.49')),
                    ('laplace_imaginary', mp.mpc(0, '0.49')),
                    ('laplace_complex', mp.mpc('0.30', '0.35'))):
                theta = scaled_theta / W
                activities = [mp.expm1(-theta * w) for w in species_weights]
                radii = [mp.expm1(abs(theta) * w) for w in species_weights]
                profiles.append((label, activities, radii, theta))
            batch_activities = [mp.expm1(a) for a in intensities]
            prefactor = mp.exp(-mp.fsum(intensities))
            for label, activities, radii, theta in profiles:
                assert all(abs(activity) <= radius_i + tolerance
                           for activity, radius_i in zip(activities, radii))
                for batch in batches:
                    assert mp.fsum(mp.asin(radii[i]) for i in batch) < mp.pi / 2
                    totals['batch_sector_checks'] += 1
                partition = graph.partition_function(activities + batch_activities)
                polynomial = absence_polynomial(activities, batch_masks, intensities)
                hard_core = prefactor * partition(graph.full)
                error = relative_error(polynomial, hard_core)
                assert error < tolerance and abs(polynomial) > tolerance
                maxima['partition_identity_relative_error'] = max(
                    maxima['partition_identity_relative_error'], error)
                totals['complex_partition_identity_checks'] += 1
                expectation = mp.fsum(p * mp.fprod(1 + activities[i] for i in range(d)
                                                   if not seen & (1 << i))
                                      for seen, p in seen_law.items())
                error = relative_error(polynomial, expectation)
                assert error < tolerance
                maxima['absence_expectation_relative_error'] = max(
                    maxima['absence_expectation_relative_error'], error)
                totals['absence_expectation_identity_checks'] += 1
                if theta is not None:
                    direct_mgf = mp.fsum(p * mp.exp(theta * scalar(x)) for x, p in pmf.items())
                    transformed = mp.exp(theta * mp.fsum(species_weights)) * polynomial
                    error = relative_error(direct_mgf, transformed)
                    assert error < tolerance
                    maxima['weighted_laplace_relative_error'] = max(
                        maxima['weighted_laplace_relative_error'], error)
                    totals['weighted_laplace_identity_checks'] += 1
                if label in ('polydisk_phases', 'laplace_complex'):
                    check_ratios(graph, partition, activities + batch_activities, radii, totals, maxima)

    gaussian_integral = mp.quad(lambda s: s * s * mp.exp(-s * s / 4), [-mp.inf, 0, mp.inf])
    assert relative_error(gaussian_integral, 4 * mp.sqrt(mp.pi)) < tolerance
    esseen_coefficient = 36 / mp.sqrt(mp.pi) + 864 / (mp.pi * mp.sqrt(2 * mp.pi))
    assert abs(esseen_coefficient - mp.mpf('130.0278279')) < mp.mpf('1e-8')
    assert esseen_coefficient < 131 < 150
    out = {
        'description': 'Deterministic supplementary checks from exhaustive finite-state sums, '
                       'independent-set sums, and fixed complex activities. These numerical checks '
                       'are not proof or formal verification of the theorem.',
        'metadata': {'python_version': platform.python_version(),
                     'python_implementation': platform.python_implementation(),
                     'platform': platform.platform(), 'mpmath_version': version('mpmath'),
                     'mpmath_decimal_precision': PRECISION, 'relative_or_absolute_tolerance': TOLERANCE,
                     'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                     'model_parameter_seed': None, 'model_selection': 'Explicit fixed parameters; no random generator.',
                     'observation_times': list(TIMES), 'maximum_cumulant_order': MAX_CUMULANT_ORDER},
        'totals': totals, 'maximum_errors_and_ratios': {key: float(value) for key, value in maxima.items()},
        'esseen_coefficient': {'formula': '36/sqrt(pi) + 864/(pi*sqrt(2*pi))',
                               'value': mp.nstr(esseen_coefficient, 40),
                               'less_than_131': True, 'less_than_150': True},
        'models': specifications, 'finite_distributions': records,
    }
    dest = args.output.expanduser().resolve()
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(json.dumps(out, indent=2) + '\n')
    print(json.dumps(totals, indent=2))
    print('All optimal-rate checks passed. Report:', dest)


if __name__ == '__main__':
    main()
