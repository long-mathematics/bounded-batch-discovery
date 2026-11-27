#!/usr/bin/env python3
"""Deterministic supplementary checks for the fixed-sample analytic results.

Requires mpmath. Uses exact rational union-state distributions, independently
enumerated independence polynomials on all induced subgraphs of small clique
incidence graphs, direct SVDs of finite Hoeffding-kernel flattenings, and exact
integer two-community occupancy counts. Complex evaluations, SVDs, and
inequality grids use 90 decimal digits. They are not interval
certificates, do not prove zero-freeness, and do not prove the conjectured
optimal fixed-sample Berry--Esseen rate. No Monte Carlo probabilities are used.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from fractions import Fraction
import hashlib
from importlib.metadata import version
from itertools import product
import json
from math import comb
from pathlib import Path
import platform

import mpmath as mp


if not __debug__:
    raise SystemExit('Verification requires assertions; run Python without -O or -OO.')

PRECISION = 90
TOLERANCE = '1e-65'
COMMUNITIES = (3, 4, 5)
COMMUNITY_SIZE_MULTIPLIERS = (1, 2, 4, 8, 16, 32, 64, 128)
TRUNCATIONS = (2, 4, 6, 8, 10)


def decimal(value, digits=45):
    return mp.nstr(value, digits)


def scalar(value):
    value = Fraction(value)
    return mp.mpf(value.numerator) / value.denominator


def weight(value):
    if value == 'sqrt(2)/4':
        return mp.sqrt(2) / 4
    if value == 'pi/8':
        return mp.pi / 8
    return scalar(value)


def members(mask):
    while mask:
        bit = mask & -mask
        yield bit.bit_length() - 1
        mask ^= bit


def batch_mask(batch):
    return sum(1 << i for i in batch)


def law(batches, raw_probabilities):
    assert len(batches) == len(raw_probabilities)
    total = sum(map(Fraction, raw_probabilities), Fraction(0))
    assert total > 0
    result = []
    for batch, probability in zip(batches, raw_probabilities):
        assert len(batch) == len(set(batch))
        probability = Fraction(probability) / total
        assert probability >= 0
        result.append((batch_mask(batch), probability))
    return result


def relative_error(left, right):
    return abs(left - right) / max(mp.mpf(1), abs(left), abs(right))


class Checks:
    def __init__(self):
        self.counts = defaultdict(int)
        self.maxima = defaultdict(lambda: mp.mpf(0))
        self.tolerance = mp.mpf(TOLERANCE)

    def identity(self, label, left, right, context=None):
        error = relative_error(left, right)
        assert error < self.tolerance, (label, context, decimal(error))
        self.counts[label] += 1
        key = label + '_relative_or_absolute_error'
        self.maxima[key] = max(self.maxima[key], error)

    def upper(self, label, value, bound, context=None):
        allowance = self.tolerance * max(mp.mpf(1), abs(value), abs(bound))
        assert value <= bound + allowance, (label, context, decimal(value), decimal(bound))
        self.counts[label] += 1
        if bound > 0:
            key = label + '_ratio'
            self.maxima[key] = max(self.maxima[key], value / bound)

    def positive(self, label, value, context=None):
        assert value > self.tolerance, (label, context, decimal(value))
        self.counts[label] += 1

    def summary(self):
        return {
            'counts': dict(sorted(self.counts.items())),
            'maximum_errors_and_ratios': {
                key: decimal(value) for key, value in sorted(self.maxima.items())},
        }


def union_distribution(draws):
    """Exact union-state dynamic programming, including nonidentical draws."""
    distribution = {0: Fraction(1)}
    for draw in draws:
        assert sum((p for _, p in draw), Fraction(0)) == 1
        updated = defaultdict(Fraction)
        for seen, probability in distribution.items():
            for batch, next_probability in draw:
                updated[seen | batch] += probability * next_probability
        distribution = {seen: p for seen, p in updated.items() if p}
    assert all(p >= 0 for p in distribution.values())
    assert sum(distribution.values(), Fraction(0)) == 1
    return distribution


def scores(weights):
    return [mp.fsum(weights[i] for i in members(mask))
            for mask in range(1 << len(weights))]


def batch_bound(draws, weights):
    return max((mp.fsum(abs(weights[i]) for i in members(batch))
                for draw in draws for batch, p in draw if p), default=mp.mpf(0))


def absence_transform(activities, draws):
    """Inclusion--exclusion expansion, independent of union-state DP."""
    total = mp.mpc(0)
    for absent_subset in range(1 << len(activities)):
        probability = Fraction(1)
        for draw in draws:
            probability *= sum((p for batch, p in draw
                                if not batch & absent_subset), Fraction(0))
        total += scalar(probability) * mp.fprod(
            activities[i] for i in members(absent_subset))
    return total


def absence_expectation(activities, distribution):
    d = len(activities)
    return mp.fsum(scalar(p) * mp.fprod(1 + activities[i] for i in range(d)
                                       if not seen & (1 << i))
                   for seen, p in distribution.items())


def mgf(distribution, state_scores, z):
    return mp.fsum(scalar(p) * mp.exp(z * state_scores[seen])
                   for seen, p in distribution.items())


class CliqueIncidenceGraph:
    """Enumerate independent subsets, then sum by an ordinary subset transform.

    No occupation-ratio recursion is used to compute any partition function.
    Every induced-subgraph partition function is obtained from the independent
    subset list, independently of the recursive inequalities being checked.
    """

    def __init__(self, d, draws):
        self.d = d
        self.option_batches = []
        self.option_probabilities = []
        self.cliques = []
        position = d
        for draw in draws:
            clique = tuple(range(position, position + len(draw)))
            self.cliques.append(clique)
            position += len(draw)
            self.option_batches.extend(batch for batch, _ in draw)
            self.option_probabilities.extend(p for _, p in draw)
        self.n = position
        self.full = (1 << self.n) - 1
        self.neighbors = [0] * self.n
        self.vertex_clique = {}
        for clique in self.cliques:
            clique_mask = sum(1 << v for v in clique)
            for v in clique:
                self.neighbors[v] |= clique_mask ^ (1 << v)
                self.vertex_clique[v] = clique_mask
        for v, batch in enumerate(self.option_batches, d):
            for i in members(batch):
                self.neighbors[v] |= 1 << i
                self.neighbors[i] |= 1 << v
        self.independent = []
        for mask in range(1 << self.n):
            vertices = tuple(members(mask))
            if all(not (self.neighbors[v] & mask) for v in vertices):
                self.independent.append((mask, vertices))

    def partitions(self, activities):
        partition = [mp.mpc(0) for _ in range(1 << self.n)]
        for mask, vertices in self.independent:
            partition[mask] = mp.fprod(activities[v] for v in vertices)
        # After this zeta transform, partition[M] is the sum over independent
        # subsets of M, i.e. the directly enumerated induced-subgraph partition.
        for bit in range(self.n):
            for mask in range(1 << self.n):
                if mask & (1 << bit):
                    partition[mask] += partition[mask ^ (1 << bit)]
        return partition


def graph_catalog():
    return [
        ('one_categorical_draw', ['1', '1', '1'],
         [law([(0,), (1,), (2,)], [1, 2, 3])]),
        ('two_iid_weighted_draws', ['1/2', '3/4', '5/4'],
         [law([(0, 1), (1, 2)], [2, 3])] * 2),
        ('two_nonidentical_draws', ['1', '1/3', '2/3', '0'],
         [law([(0, 1), (2,)], [1, 3]),
          law([(0, 2), (1, 3)], [4, 1])]),
        ('zero_activity_and_empty_batch', ['1', '1/2'],
         [law([(), (0,), (0, 1)], [1, 0, 2])]),
        ('three_cliques', ['1/4', '3/4'],
         [law([(0,), (1,)], [1, 3]), law([(), (0, 1)], [2, 1]),
          law([(0,), (1,)], [3, 2])]),
    ]


def verify_graphs(checks):
    records = []
    for name, weight_strings, draws in graph_catalog():
        weights = [weight(w) for w in weight_strings]
        W = batch_bound(draws, weights)
        graph = CliqueIncidenceGraph(len(weights), draws)
        original_distribution = union_distribution(draws)
        radii = [mp.expm1(mp.mpf('0.24') * abs(w) / W) for w in weights]
        phases = [mp.pi * (2 * i + 1) / (len(weights) + 2)
                  for i in range(len(weights))]
        profiles = [
            ('independent_complex_phases', [r * mp.exp(1j * a)
                                             for r, a in zip(radii, phases)], radii),
        ]
        z = mp.mpc('0.12', '0.18') / W
        profiles.append(('weighted_mgf', [mp.expm1(-z * w) for w in weights],
                         [mp.expm1(abs(z) * abs(w)) for w in weights]))
        for profile, left_activities, profile_radii in profiles:
            for i in range(graph.d):
                checks.upper('left_activity_radius', abs(left_activities[i]), profile_radii[i])
            for batch in graph.option_batches:
                alpha = mp.fsum(mp.asin(profile_radii[i]) for i in members(batch))
                assert alpha < mp.pi / 4
                checks.upper('right_option_sector_budget', alpha, mp.pi / 4)
            original_absence = absence_transform(left_activities, draws)
            checks.identity('graph_original_absence_dp', original_absence,
                            absence_expectation(left_activities, original_distribution), name)
            for scale in (0, 1, 7):
                right_activities = [scale * scalar(p) for p in graph.option_probabilities]
                activities = left_activities + right_activities
                partition = graph.partitions(activities)
                for mask, value in enumerate(partition):
                    checks.positive('induced_partition_nonzero_finite_evaluation', abs(value),
                                    (name, profile, scale, mask))
                    for v in members(mask):
                        without_v = mask ^ (1 << v)
                        ratio = activities[v] * partition[without_v & ~graph.neighbors[v]] / partition[without_v]
                        checks.identity('occupation_ratio_identity', value,
                                        partition[without_v] * (1 + ratio))
                        if v < graph.d:
                            checks.upper('left_occupation_radius', abs(ratio), profile_radii[v])
                            product = activities[v]
                            current = without_v
                            for neighbor in members(graph.neighbors[v] & mask):
                                child_mask = current ^ (1 << neighbor)
                                child = (activities[neighbor]
                                         * partition[child_mask & ~graph.neighbors[neighbor]]
                                         / partition[child_mask])
                                product /= 1 + child
                                current = child_mask
                            checks.identity('left_telescoping_identity', ratio, product)
                        else:
                            checks.upper('right_occupation_real_part', -mp.re(ratio), 0)
                            clique = graph.vertex_clique[v] & mask
                            H = mask & ~clique
                            q = {}
                            for option in members(clique):
                                left_neighbors = graph.neighbors[option] & ((1 << graph.d) - 1)
                                q[option] = partition[H & ~left_neighbors] / partition[H]
                                alpha = mp.fsum(mp.asin(profile_radii[i])
                                                for i in members(left_neighbors & H))
                                checks.upper('clique_deleted_ratio_sector', abs(mp.arg(q[option])), alpha)
                            denominator = 1 + mp.fsum(activities[k] * q[k]
                                                      for k in members(clique) if k != v)
                            checks.identity('right_clique_ratio_identity', ratio,
                                            activities[v] * q[v] / denominator)
                            checks.identity('clique_partition_identity', value,
                                            partition[H] * (1 + mp.fsum(activities[k] * q[k]
                                                                       for k in members(clique))))
                relaxed_draws = [[(0, Fraction(1, scale + 1))]
                                 + [(batch, p * Fraction(scale, scale + 1))
                                    for batch, p in draw] for draw in draws]
                relaxed_absence = absence_expectation(left_activities, union_distribution(relaxed_draws))
                checks.identity('normalized_clique_partition_relaxed_dp',
                                partition[graph.full] / (scale + 1)**len(draws), relaxed_absence)
            # A deterministic finite-M convergence illustration, checked
            # against the coupling error bound rather than asserting a limit.
            convergence = []
            envelope = mp.fprod(1 + r for r in profile_radii)
            for scale in (1, 10, 100, 1000):
                activities = left_activities + [scale * scalar(p) for p in graph.option_probabilities]
                normalized = graph.partitions(activities)[graph.full] / (scale + 1)**len(draws)
                error = abs(normalized - original_absence)
                bound = 2 * envelope * (1 - (mp.mpf(scale) / (scale + 1))**len(draws))
                checks.upper('finite_fugacity_coupling_bound', error, bound)
                convergence.append({'M': scale, 'absolute_error': decimal(error),
                                    'coupling_upper_bound': decimal(bound)})
            records.append({'model': name, 'profile': profile,
                            'left_vertices': graph.d, 'right_vertices': graph.n - graph.d,
                            'right_cliques': len(draws),
                            'independent_subsets_enumerated': len(graph.independent),
                            'all_induced_subgraphs_per_scale': 1 << graph.n,
                            'right_fugacity_scales': [0, 1, 7],
                            'finite_fugacity_convergence': convergence})
    return records


def finite_catalog():
    singletons = law([(0,), (1,), (2,), (3,)], [1, 3, 6, 10])
    star = law([(0, i) for i in range(1, 5)], [7, 13, 30, 50])
    return [
        *[(f'nonuniform_singletons_n{n}', ['1'] * 4, [singletons] * n)
          for n in (1, 2, 8, 40)],
        ('weighted_star', ['3/2', '1/4', '2', '0', '3/4'], [star] * 12),
        ('many_small_weights', ['1/16', '3/16', '1/8', '0', '1/4', '3/8'],
         [law([tuple(range(6)), (0, 1), (1, 2), (2, 3), (3, 4), (4, 5)],
              [2, 1, 2, 3, 4, 5])] * 5),
        ('irrational_nonnegative_weights', ['1', 'sqrt(2)/4', 'pi/8'],
         [law([(0, 1), (1, 2), (0,), ()], [3, 2, 4, 1])] * 4),
        ('nonidentical_weighted_laws', ['1/3', '2/3', '5/4', '0'],
         [law([(0, 1), (2,)], [1, 3]), law([(0, 2), (1, 3)], [4, 1]),
          law([(), (1, 2)], [2, 3]), law([(0,), (2, 3)], [1, 1])]),
        ('signed_overlapping_weights', ['1', '-1', '1/2', '-3/2'],
         [law([(0, 1), (1, 2), (2, 3), ()], [1, 3, 2, 1])] * 8),
        ('signed_nonidentical_weights', ['3/4', '-1/4', '2', '-1'],
         [law([(0, 1), (2,)], [1, 3]), law([(0, 2), (1, 3)], [4, 1]),
          law([(), (1, 2)], [2, 3]), law([(0,), (2, 3)], [1, 1])]),
        ('certain_weighted_batch', ['0', '1/2', '3/2'],
         [law([(0, 1, 2)], [1])] * 3),
        ('empty_only', ['1', '2'], [law([()], [1])] * 3),
    ]


def verify_finite_models(checks):
    records = []
    for name, weight_strings, draws in finite_catalog():
        weights = [weight(w) for w in weight_strings]
        d = len(weights)
        state_scores = scores(weights)
        distribution = union_distribution(draws)
        W = batch_bound(draws, weights)
        mean = mp.fsum(scalar(p) * state_scores[seen] for seen, p in distribution.items())
        variance = mp.fsum(scalar(p) * (state_scores[seen] - mean)**2
                          for seen, p in distribution.items())
        absent = []
        for i in range(d):
            absent.append(mp.fprod(scalar(sum((p for batch, p in draw
                                               if not batch & (1 << i)), Fraction(0)))
                                    for draw in draws))
        closed_mean = mp.fsum(w * (1 - q) for w, q in zip(weights, absent))
        closed_variance = mp.fsum(
            weights[i] * weights[j] * (
                mp.fprod(scalar(sum((p for batch, p in draw
                                    if not batch & ((1 << i) | (1 << j))), Fraction(0)))
                         for draw in draws) - absent[i] * absent[j])
            for i in range(d) for j in range(d))
        checks.identity('fixed_mean_identity', mean, closed_mean, name)
        checks.identity('fixed_variance_identity', variance, closed_variance, name)
        points = [mp.mpf('0.24'), mp.mpf('-0.24'), mp.mpc(0, '0.24'),
                  mp.mpc('0.12', '0.18'), mp.mpc('-0.12', '0.18')]
        minimum_centered_modulus = mp.inf
        for scaled_z in points:
            z = scaled_z / W if W else scaled_z
            activities = [mp.expm1(-z * w) for w in weights]
            transformed = absence_transform(activities, draws)
            checks.identity('fixed_absence_transform_dp', transformed,
                            absence_expectation(activities, distribution), name)
            direct = mgf(distribution, state_scores, z)
            checks.identity('fixed_weighted_mgf_transform', direct,
                            mp.exp(z * mp.fsum(weights)) * transformed, name)
            centered_modulus = abs(mp.exp(-z * mean) * direct)
            checks.positive('fixed_mgf_nonzero_finite_evaluation', centered_modulus, name)
            minimum_centered_modulus = min(minimum_centered_modulus, centered_modulus)
            if W:
                assert abs(z) < 1 / (4 * W)
                for draw in draws:
                    for batch, p in draw:
                        if p:
                            angle_budget = mp.fsum(mp.asin(mp.expm1(abs(z) * abs(weights[i])))
                                                   for i in members(batch))
                            checks.upper('fixed_mgf_batch_angle_budget', angle_budget, mp.pi / 4)
        records.append({
            'name': name, 'species_weights': weight_strings, 'draws': len(draws),
            'signed_weights': any(w < 0 for w in weights),
            'identically_distributed': all(draw == draws[0] for draw in draws),
            'union_states': len(distribution), 'W': decimal(W),
            'mean': decimal(mean), 'variance': decimal(variance),
            'minimum_centered_mgf_modulus_on_five_test_points': decimal(minimum_centered_modulus),
            'evaluation_scope': 'Finite complex evaluations only; not a zero-free certificate.',
        })
    return records


def verify_hoeffding_flattenings(checks):
    """Build the actual L2 kernel matrices and compute their singular values.

    Only nonnegative weights enter this proposition. Coordinates are explicit
    batches, and square roots of their exact probabilities identify L2 with
    ordinary Euclidean space. This checks the kernel, not an implementation of
    the operator upper bound. Numerical SVD is supplementary, not certified.
    """
    catalog = [
        ('low_incidence_singletons', ['1', '1/2', '1/4', '3/4'],
         law([(0,), (1,), (2,), (3,)], [1, 2, 3, 4])),
        ('mixed_incidence_overlap', ['1/2', '1'],
         law([(0,), (1,), (0, 1)], [1, 2, 1])),
    ]
    records = []
    for name, weight_strings, draw in catalog:
        weights = [weight(w) for w in weight_strings]
        d = len(weights)
        probabilities = [scalar(p) for _, p in draw]
        p_i = [mp.fsum(probabilities[a] for a, (batch, _) in enumerate(draw)
                       if batch & (1 << i)) for i in range(d)]
        q_i = [1 - p for p in p_i]
        xi = [[mp.mpf(bool(batch & (1 << i))) - p_i[i]
               for i in range(d)] for batch, _ in draw]
        covariance = [[mp.fsum(probabilities[a] * xi[a][i] * xi[a][j]
                               for a in range(len(draw)))
                       for j in range(d)] for i in range(d)]
        W = batch_bound([draw], weights)
        low_incidence = max(p_i) <= mp.mpf(1) / 2
        for n in (2, 4, 9):
            for k in range(2, min(n, 4) + 1):
                coefficients = [(-1)**(k + 1) * weights[i] * q_i[i]**(n - k)
                                for i in range(d)]

                def kernel(indices):
                    return mp.fsum(coefficients[i] * mp.fprod(xi[a][i] for a in indices)
                                   for i in range(d))

                # Test canonicality directly by integrating the final
                # coordinate, for every choice of the other coordinates.
                for prefix in product(range(len(draw)), repeat=k - 1):
                    marginal = mp.fsum(probabilities[a] * kernel(prefix + (a,))
                                       for a in range(len(draw)))
                    checks.identity('hoeffding_kernel_canonicality', marginal, 0)
                V_k = comb(n, k) * mp.fsum(
                    weights[i] * weights[j] * q_i[i]**(n - k) * q_i[j]**(n - k)
                    * covariance[i][j]**k for i in range(d) for j in range(d))
                for r in sorted({1, k // 2}):
                    left = list(product(range(len(draw)), repeat=r))
                    right = list(product(range(len(draw)), repeat=k - r))
                    left_weights = [mp.sqrt(mp.fprod(probabilities[a] for a in indices))
                                    for indices in left]
                    right_weights = [mp.sqrt(mp.fprod(probabilities[a] for a in indices))
                                     for indices in right]
                    matrix = mp.matrix([[left_weights[a] * right_weights[b]
                                         * kernel(left[a] + right[b])
                                         for b in range(len(right))] for a in range(len(left))])
                    singular_values = mp.svd_r(matrix, compute_uv=False)
                    operator_norm = singular_values[0]
                    frobenius_squared = mp.fsum(abs(matrix[a, b])**2
                                                for a in range(matrix.rows)
                                                for b in range(matrix.cols))
                    checks.identity('hoeffding_kernel_variance_identity',
                                    comb(n, k) * frobenius_squared, V_k)
                    raw_bound = W * max(p_i[i]**(mp.mpf(k) / 2) * q_i[i]**(n - k)
                                        for i in range(d))
                    checks.upper('hoeffding_flattening_operator_bound', operator_norm, raw_bound)
                    normalized = mp.sqrt(comb(n, k)) * operator_norm
                    checks.upper('hoeffding_flattening_normalized_bound', normalized, W)
                    if low_incidence:
                        checks.upper('hoeffding_flattening_geometric_decay', normalized,
                                     W * (mp.mpf(2) / 3)**(mp.mpf(k) / 2))
                    records.append({'model': name, 'n': n, 'k': k, 'split': [r, k-r],
                                    'matrix_shape': [matrix.rows, matrix.cols],
                                    'variance_component': decimal(V_k),
                                    'operator_norm': decimal(operator_norm),
                                    'operator_upper_bound': decimal(raw_bound),
                                    'normalized_operator_norm_over_W': decimal(normalized / W),
                                    'geometric_decay_hypothesis': bool(low_incidence)})
    return records


def verify_countable_truncations(checks):
    records = []
    n = 4
    z = mp.mpc('0.1', '0.2')
    # Countable singleton law P(A={i})=2^(-i-1), i>=0. Outside a
    # truncation it becomes an empty batch, not a renormalized singleton law.
    singleton_transforms = {}
    for K in TRUNCATIONS:
        draw = [(1 << i, Fraction(1, 2**(i + 1))) for i in range(K)]
        draw.append((0, Fraction(1, 2**K)))
        distribution = union_distribution([draw] * n)
        state_scores = [mp.mpf(mask.bit_count()) for mask in range(1 << K)]
        transform = mgf(distribution, state_scores, z)
        singleton_transforms[K] = transform
        expected = mp.fsum(scalar(p) * state_scores[seen] for seen, p in distribution.items())
        closed = mp.fsum(1 - (1 - mp.mpf(2)**(-i - 1))**n for i in range(K))
        checks.identity('countable_singleton_truncation_mean', expected, closed)
        checks.identity('countable_singleton_absence_identity', transform,
                        mp.exp(z * K) * absence_transform([mp.expm1(-z)] * K, [draw] * n))
        records.append({'family': 'geometric_singletons', 'truncation': K,
                        'draws': n, 'tail_probability_per_draw': str(Fraction(1, 2**K)),
                        'union_states': len(distribution), 'mean': decimal(expected)})
    reference = singleton_transforms[TRUNCATIONS[-1]]
    for K in TRUNCATIONS[:-1]:
        # Coupling yields E|D_L-D_K|<=n*2^-K; the exponential is
        # Lipschitz on [0,n] with constant |z|exp(n|z|).
        bound = abs(z) * mp.exp(n * abs(z)) * n * mp.mpf(2)**(-K)
        checks.upper('countable_singleton_transform_cauchy_bound',
                     abs(reference - singleton_transforms[K]), bound)
    # An unbounded-cardinality finite batch: A={0,...,J-1},
    # P(J=j)=2^-j for j>=1, h_i=2^(-i-1). Every batch has weight<1.
    infinite_mean = sum(((-1)**(ell + 1) * comb(n, ell)
                         * Fraction(2**ell, 2**(ell + 1) - 1)
                         for ell in range(1, n + 1)), Fraction(0))
    prefix_transforms = {}
    for K in TRUNCATIONS:
        draw = [((1 << j) - 1, Fraction(1, 2**j)) for j in range(1, K)]
        draw.append(((1 << K) - 1, Fraction(1, 2**(K - 1))))
        distribution = union_distribution([draw] * n)
        rational_scores = [sum((Fraction(1, 2**(i + 1)) for i in members(mask)), Fraction(0))
                           for mask in range(1 << K)]
        expected = sum((p * rational_scores[seen] for seen, p in distribution.items()), Fraction(0))
        exact_tail = sum(((-1)**(ell + 1) * comb(n, ell)
                          * Fraction(2**ell, (2**(ell + 1) - 1) * 2**((ell + 1) * K))
                          for ell in range(1, n + 1)), Fraction(0))
        assert infinite_mean - expected == exact_tail >= 0
        checks.counts['countable_prefix_exact_mean_tail_identity'] += 1
        checks.upper('countable_prefix_mean_tail_bound', scalar(exact_tail),
                     scalar(Fraction(2 * n, 3 * 4**K)))
        transform = mgf(distribution, [scalar(s) for s in rational_scores], z)
        prefix_transforms[K] = transform
        activities = [mp.expm1(-z * mp.mpf(2)**(-i - 1)) for i in range(K)]
        checks.identity('countable_prefix_absence_identity', transform,
                        mp.exp(z * (1 - mp.mpf(2)**(-K)))
                        * absence_transform(activities, [draw] * n))
        records.append({'family': 'geometric_prefixes_with_unbounded_batch_cardinality',
                        'truncation': K, 'draws': n, 'union_states': len(distribution),
                        'batch_weight_bound': str(Fraction(2**K - 1, 2**K)),
                        'mean_exact': str(expected), 'infinite_model_mean_exact': str(infinite_mean),
                        'exact_mean_tail': str(exact_tail)})
    reference = prefix_transforms[TRUNCATIONS[-1]]
    for K in TRUNCATIONS[:-1]:
        bound = abs(z) * mp.exp(abs(z)) * mp.mpf(2)**(-K)
        checks.upper('countable_prefix_transform_cauchy_bound',
                     abs(reference - prefix_transforms[K]), bound)
    return records


def verify_analytic_parameters(checks):
    records = []
    # This grid illustrates the elementary analytic inequality in the proof;
    # its continuum validity follows from the manuscript, not from the grid.
    for numerator in range(251):
        a = mp.mpf(numerator) / 1000
        checks.upper('arcsin_exponential_grid', mp.asin(mp.expm1(a)), mp.pi * a)
    R = mp.exp(mp.mpf(1) / 16)
    r = mp.exp(mp.mpf(1) / 32)
    for n in (1, 2, 3, 10, 100, 10**6, 10**12, 10**100):
        delta = 1 / (64 * mp.log(2 * n + 2))
        epsilon = delta / 4
        log_radius = mp.sqrt((mp.mpf(1) / 16)**2 + delta**2)
        checks.upper('harnack_sector_inside_log_disk', log_radius, mp.mpf(1) / 8)
        attenuation = (r / R)**(1 / delta)
        checks.identity('ms_boundary_attenuation_identity', attenuation,
                        1 / mp.mpf(2 * n + 2)**2)
        boundary = attenuation * n / 2
        checks.upper('ms_boundary_smallness', boundary, mp.mpf(1) / 32)
        checks.upper('ms_boundary_lemma_threshold', boundary, mp.mpf(3) / 8)
        checks.upper('ms_disk_angular_containment', mp.asin(epsilon), delta / 2)
        checks.upper('ms_disk_inner_radial_containment', 1 / r, 1 - epsilon)
        checks.upper('ms_disk_outer_radial_containment', 1 + epsilon, r)
        sigma = 65 / epsilon
        variance_term = sigma**2 / 2 * (epsilon / 32)**2
        checks.positive('ms_corollary_variance_series_threshold', variance_term - 1)
        checks.identity('harnack_disk_factor',
                        (mp.mpf(1) / 4 + mp.mpf(1) / 8)
                        / (mp.mpf(1) / 4 - mp.mpf(1) / 8), 3)
        records.append({'n': str(n), 'delta': decimal(delta), 'epsilon': decimal(epsilon),
                        'maximum_log_modulus_on_sector': decimal(log_radius),
                        'boundary_quantity_upper_bound': decimal(boundary),
                        'variance_term_at_sigma_65_over_epsilon': decimal(variance_term)})
    return records


def choose(n, k):
    return comb(n, k) if 0 <= k <= n else 0


def community_integer_dp(M, sample_sizes):
    """Count ordered sequences of pair batches exactly, using symmetry."""
    counts = {(0, 0): 1}
    selected = {}
    for n in range(1, max(sample_sizes) + 1):
        updated = defaultdict(int)
        for (a, b), count in counts.items():
            for new in range(3):
                multiplicity_a = choose(M - a, new) * choose(a, 2 - new)
                multiplicity_b = choose(M - b, new) * choose(b, 2 - new)
                if multiplicity_a:
                    updated[a + new, b] += count * multiplicity_a
                if multiplicity_b:
                    updated[a, b + new] += count * multiplicity_b
        counts = dict(updated)
        if n in sample_sizes:
            by_degree = [0] * (2 * M + 1)
            for (a, b), count in counts.items():
                by_degree[a + b] += count
            selected[n] = by_degree
    return selected


def community_inclusion_exclusion(M, n):
    """Independent exact integer count for each specified union, then aggregate."""
    by_degree = [0] * (2 * M + 1)
    for a in range(M + 1):
        for b in range(M + 1):
            exact_union_count = sum(
                (-1)**(a - u + b - v) * comb(a, u) * comb(b, v)
                * (choose(u, 2) + choose(v, 2))**n
                for u in range(a + 1) for v in range(b + 1))
            assert exact_union_count >= 0
            by_degree[a + b] += comb(M, a) * comb(M, b) * exact_union_count
    return by_degree


def verify_communities(checks):
    records = []
    for M in COMMUNITIES:
        denominator_base = 2 * comb(M, 2)
        sample_sizes = [M * multiplier for multiplier in COMMUNITY_SIZE_MULTIPLIERS]
        exact_dp = community_integer_dp(M, sample_sizes)
        for a in range(M + 1):
            for b in range(M + 1):
                capacity = scalar(Fraction(choose(a, 2) + choose(b, 2), denominator_base))
                bound = mp.power(2, mp.mpf(a + b) / M - 2)
                checks.upper('two_community_capacity_inequality', capacity, bound)
                if (a, b) in ((M, 0), (0, M), (M, M)):
                    checks.identity('two_community_capacity_equality_cases', capacity, bound)
                else:
                    checks.positive('two_community_capacity_strict_cases', bound - capacity)
        rows = []
        for n in sample_sizes:
            counts = community_inclusion_exclusion(M, n)
            assert counts == exact_dp[n]
            checks.counts['two_community_integer_dp_inclusion_exclusion_identity'] += 1
            denominator = denominator_base**n
            assert sum(counts) == denominator
            checks.counts['two_community_exact_probability_normalization'] += 1
            scaled = [mp.power(2, mp.mpf(n) * (2 - mp.mpf(d) / M))
                      * mp.mpf(count) / denominator for d, count in enumerate(counts)]
            target = [mp.mpf(2) if d == M else mp.mpf(1) if d == 2 * M else mp.mpf(0)
                      for d in range(2 * M + 1)]
            maximum_error = max(abs(value - limit) for value, limit in zip(scaled, target))
            off_degree = max(scaled[d] for d in range(2 * M + 1) if d not in (M, 2 * M))
            rows.append({'n': n, 'coefficient_at_M': decimal(scaled[M]),
                         'coefficient_at_2M': decimal(scaled[2 * M]),
                         'maximum_other_scaled_coefficient': decimal(off_degree),
                         'maximum_absolute_coefficient_error': decimal(maximum_error),
                         'scaled_coefficients_in_increasing_degree': list(map(decimal, scaled))})
        assert mp.mpf(rows[-1]['maximum_absolute_coefficient_error']) < mp.mpf('1e-20')
        checks.counts['two_community_final_coefficient_convergence_illustration'] += 1
        records.append({'M': M, 'B': 2, 'batch_options': denominator_base,
                        'rescaled_polynomial': '2^(2n) G_(M,n)(2^(-n/M) z)',
                        'proved_limit_in_manuscript': f'z^{2*M} + 2 z^{M}',
                        'arithmetic': 'Exact integer DP and independent inclusion--exclusion; decimal rescaling only.',
                        'interpretation': 'Finite illustrations of the proved coefficient limit; no numerical root certificate.',
                        'cases': rows})
    return records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path,
                        default=Path(__file__).resolve().parents[1] / 'reports' / 'fixed_sample_verification.json')
    args = parser.parse_args()
    mp.mp.dps = PRECISION
    checks = Checks()
    graphs = verify_graphs(checks)
    finite = verify_finite_models(checks)
    flattenings = verify_hoeffding_flattenings(checks)
    countable = verify_countable_truncations(checks)
    parameters = verify_analytic_parameters(checks)
    communities = verify_communities(checks)
    script = Path(__file__).resolve()
    report = {
        'status': 'passed',
        'description': 'Deterministic supplementary checks for the fixed-sample zero-free disk, '
                       'logarithmic-rate analytic parameters, and global-sector obstruction. '
                       'The manuscript proofs establish the stated results. These checks are '
                       'not formal or interval-certified verification and do not prove the '
                       'optimal fixed-sample Berry--Esseen conjecture.',
        'metadata': {
            'python_version': platform.python_version(),
            'python_implementation': platform.python_implementation(),
            'platform': platform.platform(), 'mpmath_version': version('mpmath'),
            'mpmath_decimal_precision': PRECISION,
            'high_precision_relative_or_absolute_tolerance': TOLERANCE,
            'script_sha256': hashlib.sha256(script.read_bytes()).hexdigest(),
            'manuscript_source_path': 'papers/fixed_sample_bounds.tex',
            'model_parameter_seed': None,
            'model_selection': 'Explicit rational laws and specified real weights; no random generator.',
            'fixed_sample_arithmetic': 'Exact fractions for union-state probabilities; real weighted transforms at 90 digits.',
            'graph_arithmetic': 'Independent-subset enumeration followed by subset summation; all small induced subgraphs checked.',
            'hoeffding_arithmetic': 'Explicit probability-weighted kernel matrices, canonicality sums, and mpmath SVD; nonnegative weights only.',
            'two_community_arithmetic': 'Exact integer sequence counts from two independent methods, then 90-digit rescaling.',
            'zero_free_scope': 'Finite evaluations are non-certificates; continuum zero-freeness requires the proof.',
        },
        'case_totals': {
            'clique_graphs': len(graph_catalog()), 'clique_graph_profiles': len(graphs),
            'finite_distribution_models': len(finite),
            'nonidentical_finite_models': sum(not row['identically_distributed'] for row in finite),
            'signed_weight_finite_models': sum(row['signed_weights'] for row in finite),
            'hoeffding_flattening_svd_cases': len(flattenings),
            'countable_finite_truncations': len(countable),
            'analytic_parameter_cases': len(parameters),
            'two_community_models': len(communities),
            'two_community_sample_size_cases': sum(len(row['cases']) for row in communities),
        },
        **checks.summary(),
        'clique_graph_checks': graphs, 'finite_distribution_checks': finite,
        'hoeffding_flattening_svd_checks': flattenings,
        'countable_truncation_checks': countable,
        'harnack_and_michelen_sahasrabudhe_parameters': parameters,
        'two_community_coefficient_checks': communities,
    }
    destination = args.output.expanduser().resolve()
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(dict(sorted(checks.counts.items())), indent=2))
    print('All fixed-sample checks passed. Report:', destination)


if __name__ == '__main__':
    main()
