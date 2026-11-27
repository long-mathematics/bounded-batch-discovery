#!/usr/bin/env python3
"""Deterministic supplementary checks for quantitative discovery extensions.

Requires mpmath and scipy. Uses 80 decimal digits, exhaustive finite Poisson
presence states, exact rational fixed-sample union laws, independent quadrature,
and finite-binomial cross-checks. Large-binomial examples use a scaled beta
integral; the analytic bound on the omitted integral is reported separately
from the numerical quadrature error. These checks supplement the proofs and
are neither formal verification nor certified interval computations.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from fractions import Fraction
import hashlib
from importlib.metadata import version
import json
from math import comb
from pathlib import Path
import platform

import mpmath as mp
from scipy.stats import binom

from verify_optimal_rate import (
    absence_polynomial,
    members,
    models as original_models,
    presence_distribution,
    relative_error,
    scalar,
)


if not __debug__:
    raise SystemExit('Verification requires assertions; run Python without -O or -OO.')

PRECISION = 80
TOLERANCE = '1e-60'
DOUBLE_TOLERANCE = '1e-10'
TIMES = ('0.05', '1', '8')
SAMPLE_SIZES = (1, 2, 3, 10, 40)
TILT_PROFILES = ('0', '1/4096', '1/2048', '1/16')
ENVELOPE_THRESHOLDS = (1, 2, 3, 5, 10, 25)
BETA_CUTOFF = 256


def decimal(value, digits=40):
    return mp.nstr(value, digits)


def model_catalog():
    return original_models() + [
        ('certain_singleton', [(0,)], [1], ['1']),
        ('empty_only', [()], [1], ['1', '2']),
        ('certain_weighted_batch', [(0, 1, 2)], [1], ['0', '1/2', '3/2']),
        ('two_uniform_singletons', [(0,), (1,)], [1, 1], ['1', '1']),
        ('rare_weighted_overlap', [(0,), (1,), (1, 2)], [1, 499, 500],
         ['1', '3/10', '17/10']),
    ]


def prepare_model(specification):
    name, batches, raw_rates, weights = specification
    assert len(batches) == len(raw_rates) and sum(raw_rates) > 0
    assert all(rate >= 0 for rate in raw_rates)
    assert all(len(set(batch)) == len(batch) for batch in batches)
    probabilities = [Fraction(rate, sum(raw_rates)) for rate in raw_rates]
    rational_weights = [Fraction(w) for w in weights]
    assert all(w >= 0 for w in rational_weights)
    d = len(weights)
    assert all(all(0 <= i < d for i in batch) for batch in batches)
    batch_masks = [sum(1 << i for i in batch) for batch in batches]
    scores = [sum((rational_weights[i] for i in members(mask)), Fraction(0))
              for mask in range(1 << d)]
    rates = [sum((p for mask, p in zip(batch_masks, probabilities)
                  if mask & (1 << i)), Fraction(0)) for i in range(d)]
    pair_rates = [[sum((p for mask, p in zip(batch_masks, probabilities)
                       if mask & (1 << i) and mask & (1 << j)), Fraction(0))
                   for j in range(d)] for i in range(d)]
    positive_masks = [mask for mask, p in zip(batch_masks, probabilities) if p > 0]
    W = max((scores[mask] for mask in positive_masks), default=Fraction(0))
    return {
        'name': name, 'batches': batches, 'raw_rates': raw_rates,
        'probabilities': probabilities, 'weight_strings': weights,
        'weights': rational_weights, 'd': d, 'batch_masks': batch_masks,
        'scores': scores, 'rates': rates, 'pair_rates': pair_rates, 'W': W,
        'mu1': sum((p * scores[mask] for mask, p in zip(batch_masks, probabilities)),
                   Fraction(0)),
        'mu2': sum((p * scores[mask]**2 for mask, p in zip(batch_masks, probabilities)),
                   Fraction(0)),
    }


class Checks:
    def __init__(self):
        self.counts = defaultdict(int)
        self.maxima = defaultdict(lambda: mp.mpf(0))
        self.tolerance = mp.mpf(TOLERANCE)

    def identity(self, label, left, right, context=None):
        error = relative_error(left, right)
        assert error < self.tolerance, (label, context, decimal(error))
        self.counts[label] += 1
        self.maxima[label + '_relative_or_absolute_error'] = max(
            self.maxima[label + '_relative_or_absolute_error'], error)

    def upper(self, label, value, bound, context=None):
        allowance = self.tolerance * max(mp.mpf(1), abs(value), abs(bound))
        assert value <= bound + allowance, (label, context, decimal(value), decimal(bound))
        self.counts[label] += 1
        if bound > 0:
            self.maxima[label + '_ratio'] = max(self.maxima[label + '_ratio'], value / bound)

    def summary(self):
        return {
            'counts': dict(sorted(self.counts.items())),
            'maximum_errors_and_ratios': {key: decimal(value)
                                          for key, value in sorted(self.maxima.items())},
        }


def score_distribution(model, union_law):
    pmf = defaultdict(lambda: mp.mpf(0))
    for mask, probability in union_law.items():
        pmf[model['scores'][mask]] += probability
    return {score: probability for score, probability in sorted(pmf.items()) if probability > 0}


def finite_moments(pmf):
    mean = mp.fsum(scalar(score) * probability for score, probability in pmf.items())
    variance = mp.fsum((scalar(score) - mean)**2 * probability
                       for score, probability in pmf.items())
    return mean, variance


def normal_cdf(x):
    return mp.erfc(-x / mp.sqrt(2)) / 2


def finite_kolmogorov(points):
    cumulative = mp.mpf(0)
    distance = mp.mpf(0)
    for value, probability in sorted(points):
        gaussian = normal_cdf(value)
        distance = max(distance, abs(cumulative - gaussian),
                       abs(cumulative + probability - gaussian))
        cumulative += probability
    return distance


def verify_tilts(catalog, checks):
    records = []
    # The nine original models are deliberately reused so that the new checks
    # are directly comparable with the earlier optimal-rate report.
    for model in catalog[:len(original_models())]:
        W = scalar(model['W'])
        assert W > 0
        weights = [scalar(w) for w in model['weights']]
        for time_string in TIMES:
            t = scalar(time_string)
            intensities = [t * scalar(p) for p in model['probabilities']]
            law = presence_distribution(model['batch_masks'], intensities)
            pmf = score_distribution(model, law)
            mean, variance = finite_moments(pmf)
            sigma = mp.sqrt(variance)
            delta = W / sigma
            checks.identity('poisson_distribution_normalization', mp.fsum(pmf.values()), 1)

            def transformed_mgf(theta):
                activities = [mp.expm1(-theta * w) for w in weights]
                return (mp.exp(theta * mp.fsum(weights))
                        * absence_polynomial(activities, model['batch_masks'], intensities))

            for sign in (-1, 1):
                points = [(sign * (scalar(score) - mean) / sigma, probability)
                          for score, probability in pmf.items()]
                for profile in TILT_PROFILES:
                    # h*delta is fixed. The physical tilt of X is sign*profile/W.
                    scaled_h = scalar(profile)
                    h = scaled_h / delta
                    theta = sign * scaled_h / W
                    direct_laplace = mp.fsum(p * mp.exp(theta * scalar(score))
                                             for score, p in pmf.items())
                    checks.identity('tilted_laplace_identity', direct_laplace,
                                    transformed_mgf(theta), (model['name'], time_string, sign, profile))
                    tilted_scores = [(scalar(score), p * mp.exp(theta * scalar(score)) / direct_laplace)
                                     for score, p in pmf.items()]
                    tilted_mean = mp.fsum(y * p for y, p in tilted_scores)
                    tilted_variance = mp.fsum((y - tilted_mean)**2 * p for y, p in tilted_scores)
                    log_transform = lambda z: mp.log(transformed_mgf(z))
                    checks.identity('tilted_mean_derivative_identity', tilted_mean,
                                    mp.diff(log_transform, theta), model['name'])
                    checks.identity('tilted_variance_derivative_identity', tilted_variance,
                                    mp.diff(log_transform, theta, 2), model['name'])

                    K = mp.log(mp.fsum(p * mp.exp(h * y) for y, p in points))
                    tilted = [(y, p * mp.exp(h * y - K)) for y, p in points]
                    checks.identity('tilted_distribution_normalization', mp.fsum(p for _, p in tilted), 1)
                    x = mp.fsum(y * p for y, p in tilted)
                    v = mp.fsum((y - x)**2 * p for y, p in tilted)
                    checks.identity('standardized_tilt_mean', x, sign * (tilted_mean - mean) / sigma)
                    checks.identity('standardized_tilt_variance', v, tilted_variance / variance)
                    a = h * mp.sqrt(v)
                    standardized = [((y - x) / mp.sqrt(v), p) for y, p in tilted]
                    distance = finite_kolmogorov(standardized)
                    in_tilt_domain = scaled_h <= mp.mpf(1) / 1024
                    if in_tilt_domain:
                        checks.upper('tilted_variance_lower_bound', mp.mpf(19) / 20, v)
                        checks.upper('tilted_variance_upper_bound', v, mp.mpf(21) / 20)
                        checks.upper('tilted_kolmogorov_bound', distance, 200 * delta)

                    for strict in (False, True):
                        condition = (lambda y, threshold: y > threshold) if strict else (
                            lambda y, threshold: y >= threshold)
                        direct_tail = mp.fsum(p for y, p in points if condition(y, x))
                        laplace_tail = mp.fsum(p * mp.exp(-a * ((y - x) / mp.sqrt(v)))
                                              for y, p in tilted if condition(y, x))
                        checks.identity('undo_tilt_tail_identity', direct_tail,
                                        mp.exp(K - h * x) * laplace_tail)

                        # An actual support point checks endpoint atoms, even
                        # when it is not the tilted mean required by the theorem.
                        endpoint = sorted(y for y, _ in points)[len(points) // 2]
                        direct_endpoint_tail = mp.fsum(p for y, p in points if condition(y, endpoint))
                        tilted_endpoint_tail = mp.fsum(p * mp.exp(-h * (y - endpoint))
                                                       for y, p in tilted if condition(y, endpoint))
                        checks.identity('undo_tilt_support_endpoint_identity', direct_endpoint_tail,
                                        mp.exp(K - h * endpoint) * tilted_endpoint_tail)
                    direct_atom = mp.fsum(p for y, p in points if y == endpoint)
                    tilted_atom = mp.fsum(p for y, p in tilted if y == endpoint)
                    checks.identity('undo_tilt_atom_identity', direct_atom,
                                    mp.exp(K - h * endpoint) * tilted_atom)
                    records.append({
                        'model': model['name'], 'time': time_string, 'sign': sign,
                        'h_times_delta': profile, 'delta': decimal(delta),
                        'tilted_standardized_mean': decimal(x),
                        'tilted_standardized_variance': decimal(v),
                        'tilted_kolmogorov_distance': decimal(distance),
                        'inside_proved_tilt_domain': in_tilt_domain,
                        'kolmogorov_bound_nonvacuous': bool(200 * delta < 1),
                    })
    return records


def exact_fixed_union_laws(model):
    law = {0: Fraction(1)}
    snapshots = {}
    for n in range(1, max(SAMPLE_SIZES) + 1):
        next_law = defaultdict(Fraction)
        for union, probability in law.items():
            for mask, p in zip(model['batch_masks'], model['probabilities']):
                next_law[union | mask] += probability * p
        law = {mask: p for mask, p in next_law.items() if p}
        assert sum(law.values(), Fraction(0)) == 1
        if n in SAMPLE_SIZES:
            snapshots[n] = law
    return snapshots


def joint_poisson_moments(model, n):
    """Enumerate presence states and conditional Poisson count moments.

    This supplies E[X N] and E[(X-fN)^2] independently of the pair covariance
    formulas. Empty batches contribute to N and not to X.
    """
    intensities = [n * scalar(p) for p in model['probabilities']]
    absent = [mp.exp(-a) for a in intensities]
    present = [-mp.expm1(-a) for a in intensities]
    count_means, count_variances = [], []
    for a, p in zip(intensities, present):
        if a == 0:
            count_means.append(mp.mpf(0))
            count_variances.append(mp.mpf(0))
        else:
            conditional_mean = a / p
            count_means.append(conditional_mean)
            count_variances.append((a + a*a) / p - conditional_mean**2)
    result = {key: mp.mpf(0) for key in ('mass', 'X', 'X2', 'N', 'N2', 'XN')}
    for state in range(1 << len(model['batch_masks'])):
        probability = mp.mpf(1)
        union = 0
        mean_n = mp.mpf(0)
        variance_n = mp.mpf(0)
        for j, mask in enumerate(model['batch_masks']):
            if state & (1 << j):
                probability *= present[j]
                union |= mask
                mean_n += count_means[j]
                variance_n += count_variances[j]
            else:
                probability *= absent[j]
        score = scalar(model['scores'][union])
        result['mass'] += probability
        result['X'] += probability * score
        result['X2'] += probability * score**2
        result['N'] += probability * mean_n
        result['N2'] += probability * (variance_n + mean_n**2)
        result['XN'] += probability * score * mean_n
    return result


def verify_fixed_samples(catalog, checks):
    records = []
    for model in catalog:
        fixed_laws = exact_fixed_union_laws(model)
        d = model['d']
        weights = [scalar(w) for w in model['weights']]
        rates = [scalar(rate) for rate in model['rates']]
        pair = [[scalar(c) for c in row] for row in model['pair_rates']]
        mu1, mu2 = scalar(model['mu1']), scalar(model['mu2'])
        for n in SAMPLE_SIZES:
            fixed = fixed_laws[n]
            exact_mean = sum((p * model['scores'][mask] for mask, p in fixed.items()), Fraction(0))
            exact_variance = sum((p * (model['scores'][mask] - exact_mean)**2
                                  for mask, p in fixed.items()), Fraction(0))
            fixed_mean, fixed_variance = scalar(exact_mean), scalar(exact_variance)
            intensities = [n * scalar(p) for p in model['probabilities']]
            poisson_law = presence_distribution(model['batch_masks'], intensities)
            poisson_mean, poisson_variance = finite_moments(score_distribution(model, poisson_law))
            context = (model['name'], n)
            closed_fixed_mean = mp.fsum(weights[i] * (1 - (1 - rates[i])**n) for i in range(d))
            closed_fixed_variance = mp.fsum(
                weights[i] * weights[j]
                * ((1 - rates[i] - rates[j] + pair[i][j])**n
                   - (1 - rates[i])**n * (1 - rates[j])**n)
                for i in range(d) for j in range(d))
            closed_poisson_variance = mp.fsum(
                weights[i] * weights[j] * mp.exp(-n * (rates[i] + rates[j]))
                * mp.expm1(n * pair[i][j]) for i in range(d) for j in range(d))
            checks.identity('fixed_mean_covariance_formula', fixed_mean, closed_fixed_mean, context)
            checks.identity('fixed_variance_covariance_formula', fixed_variance, closed_fixed_variance, context)
            checks.identity('poisson_variance_covariance_formula', poisson_variance,
                            closed_poisson_variance, context)
            f = mp.fsum(weights[i] * rates[i] * mp.exp(-n * rates[i]) for i in range(d))
            tau2 = poisson_variance - n * f*f
            G = [mp.fsum(weights[i] * mp.exp(-n * rates[i]) for i in batch)
                 for batch in model['batches']]
            G_mean = mp.fsum(scalar(p) * value for p, value in zip(model['probabilities'], G))
            positive = n * mp.fsum(scalar(p) * (value - G_mean)**2
                                    for p, value in zip(model['probabilities'], G))
            positive += mp.fsum(
                weights[i] * weights[j] * mp.exp(-n * (rates[i] + rates[j]))
                * (mp.expm1(n * pair[i][j]) - n * pair[i][j])
                for i in range(d) for j in range(d))
            checks.identity('positive_residual_variance_identity', tau2, positive, context)
            checks.upper('residual_variance_nonnegative', -tau2, 0, context)

            joint = joint_poisson_moments(model, n)
            checks.identity('joint_poisson_normalization', joint['mass'], 1, context)
            checks.identity('joint_poisson_first_count_moment', joint['N'], n, context)
            checks.identity('joint_poisson_second_count_moment', joint['N2'], n + n*n, context)
            checks.identity('joint_poisson_score_mean', joint['X'], poisson_mean, context)
            covariance = joint['XN'] - joint['X'] * joint['N']
            checks.identity('poisson_projection_covariance', covariance, n * f, context)
            projected_mean = joint['X'] - f * joint['N']
            projected_variance = (joint['X2'] - 2*f*joint['XN'] + f*f*joint['N2']
                                  - projected_mean**2)
            checks.identity('poisson_projection_variance', projected_variance, tau2, context)

            mean_difference = fixed_mean - poisson_mean
            remainder = fixed_variance - tau2
            upper = 2 * (mu2 + mu1*mu1)
            lower_magnitude = upper + 2 * mp.exp(-2) * mu1*mu1
            checks.upper('fixed_mean_difference_nonnegative', -mean_difference, 0, context)
            checks.upper('fixed_mean_difference_upper_bound', mean_difference, 2 * mu1, context)
            checks.upper('fixed_variance_correction_upper_bound', remainder, upper, context)
            checks.upper('fixed_variance_correction_lower_bound', -remainder, lower_magnitude, context)
            checks.upper('fixed_variance_absolute_correction', abs(remainder),
                         2*mu2 + (2 + 2*mp.exp(-2))*mu1*mu1, context)
            records.append({
                'model': model['name'], 'sample_size': n, 'union_support_size': len(fixed),
                'fixed_mean': decimal(fixed_mean), 'fixed_variance': decimal(fixed_variance),
                'poisson_mean': decimal(poisson_mean), 'poisson_variance': decimal(poisson_variance),
                'f': decimal(f), 'residual_variance': decimal(tau2),
                'mean_difference': decimal(mean_difference),
                'variance_correction_remainder': decimal(remainder),
                'absolute_remainder_bound': decimal(lower_magnitude),
            })
    return records


def poisson_mass(x, k):
    if x == 0:
        return mp.mpf(1) if k == 0 else mp.mpf(0)
    return mp.exp(-x + k * mp.log(x) - mp.loggamma(k + 1))


def poisson_envelope(x, r):
    return max(poisson_mass(x, k) for k in range(r))


def ell(r):
    return mp.mpf(1) + mp.fsum(poisson_mass(mp.mpf(k), k) for k in range(1, r))


def envelope_integral_gamma(r):
    return (mp.fsum(mp.gammainc(k + 1, k, k + 1) / mp.factorial(k)
                    for k in range(r - 1))
            + mp.gammainc(r, r - 1, mp.inf) / mp.factorial(r - 1))


def two_rate_integral(r, a, b):
    """Integrate the piecewise polynomial-exponential product exactly by gamma functions."""
    a, b = Fraction(a), Fraction(b)
    endpoints = sorted({Fraction(0)} | {Fraction(k) / a for k in range(1, r)}
                       | {Fraction(k) / b for k in range(1, r)})
    rate_a, rate_b = scalar(a), scalar(b)
    values = []
    for j, left in enumerate(endpoints):
        right = endpoints[j + 1] if j + 1 < len(endpoints) else None
        if right is None:
            ka = kb = r - 1
        else:
            midpoint = (left + right) / 2
            ka, kb = min(int(a * midpoint), r - 1), min(int(b * midpoint), r - 1)
        power = ka + kb
        coefficient = (rate_a**ka * rate_b**kb
                       / (mp.factorial(ka) * mp.factorial(kb) * (rate_a + rate_b)**(power + 1)))
        values.append(coefficient * mp.gammainc(
            power + 1, (rate_a + rate_b) * scalar(left),
            mp.inf if right is None else (rate_a + rate_b) * scalar(right)))
    return mp.fsum(values), endpoints


def verify_envelopes_and_loads(catalog, checks):
    envelope_records, product_records, load_records = [], [], []
    for r in ENVELOPE_THRESHOLDS:
        value = ell(r)
        integral = envelope_integral_gamma(r)
        checks.identity('poisson_envelope_gamma_integral', integral, value, r)
        upper = 1 + mp.sqrt(2 * mp.mpf(r - 1) / mp.pi)
        checks.upper('poisson_envelope_explicit_upper_bound', value, upper, r)
        for k in range(r + 2):
            for offset in (Fraction(0), Fraction(1, 4), Fraction(3, 4)):
                x = scalar(Fraction(k) + offset)
                mode = min(int(mp.floor(x)), r - 1)
                checks.identity('poisson_envelope_mode_identity', poisson_envelope(x, r),
                                poisson_mass(x, mode), (r, str(x)))
        quadrature = None
        if r <= 10:
            quadrature = mp.fsum(
                mp.quad(lambda x, rr=r: poisson_envelope(x, rr), [k, k + 1])
                for k in range(r - 1))
            quadrature += mp.quad(lambda x: poisson_envelope(x, r), [r - 1, mp.inf])
            checks.identity('poisson_envelope_direct_quadrature', quadrature, value, r)
        envelope_records.append({'r': r, 'ell_r': decimal(value),
                                 'explicit_upper_bound': decimal(upper),
                                 'piecewise_gamma_integral': decimal(integral),
                                 'direct_quadrature': None if quadrature is None else decimal(quadrature)})

        previous = mp.mpf(0)
        for a, b in [('1', '1'), ('1', '1/2'), ('1', '1/10'), ('1', '1/10000'), ('2/3', '7/5')]:
            product, endpoints = two_rate_integral(r, a, b)
            bound = value / max(scalar(a), scalar(b))
            checks.upper('two_rate_envelope_bound', product, bound, (r, a, b))
            if a == '1':
                checks.upper('two_rate_small_rate_monotonicity', previous, product, (r, a, b))
                previous = product
            direct = None
            if r <= 3 and (a, b) in [('1', '1/2'), ('2/3', '7/5')]:
                intervals = [scalar(x) for x in endpoints] + [mp.inf]
                direct = mp.quad(lambda u: poisson_envelope(scalar(a) * u, r)
                                * poisson_envelope(scalar(b) * u, r), intervals)
                checks.identity('two_rate_envelope_direct_quadrature', direct, product, (r, a, b))
            if r == 1:
                checks.identity('threshold_one_product_identity', product,
                                1 / (scalar(a) + scalar(b)))
            product_records.append({'r': r, 'a': a, 'b': b, 'integral': decimal(product),
                                    'upper_bound': decimal(bound), 'ratio_to_bound': decimal(product / bound),
                                    'direct_quadrature': None if direct is None else decimal(direct)})

    for model in catalog:
        W = scalar(model['W'])
        if W == 0:
            continue
        weights = [scalar(w) for w in model['weights']]
        rates = [scalar(rate) for rate in model['rates']]
        pair = [[scalar(c) for c in row] for row in model['pair_rates']]
        for i in range(model['d']):
            row = mp.fsum(weights[j] * pair[i][j] for j in range(model['d']))
            checks.upper('weighted_incidence_row_bound', row, W * rates[i], (model['name'], i))
        for r in (1, 2, 5, 25):
            L = ell(r)
            coefficients = [[(L * weights[i] * weights[j] * pair[i][j]
                              / max(rates[i], rates[j]) if max(rates[i], rates[j]) > 0 else mp.mpf(0))
                             for j in range(model['d'])] for i in range(model['d'])]
            for batch, probability in zip(model['batches'], model['probabilities']):
                if probability == 0:
                    continue
                load = mp.fsum(coefficients[i][j] for i in range(model['d']) for j in range(model['d'])
                               if i in batch or j in batch)
                bound = 2 * L * W*W
                checks.upper('threshold_weighted_coordinate_load', load, bound, (model['name'], r, batch))
                load_records.append({'model': model['name'], 'r': r, 'batch': batch,
                                     'coefficient_load': decimal(load), 'upper_bound': decimal(bound),
                                     'ratio_to_bound': decimal(load / bound)})
    asymptotic_records = [
        {'r': r, 'ell_r': decimal(ell(r)),
         'ell_r_over_sqrt_2r_over_pi': decimal(ell(r) / mp.sqrt(2 * mp.mpf(r) / mp.pi)),
         'interpretation': 'Finite illustration of the proved asymptotic, not a convergence test.'}
        for r in (10, 100, 1000)]
    return {'envelope_integrals': envelope_records, 'two_rate_integrals': product_records,
            'weighted_loads': load_records, 'asymptotic_illustrations': asymptotic_records}


def verify_constants(checks):
    parameters = {
        'tilted_esseen_constant': 52 / mp.sqrt(mp.pi) + 1248 / (mp.pi * mp.sqrt(2 * mp.pi)),
        'tilt_path_radius_times_delta': mp.mpf(1) / 1024 + mp.sqrt(mp.mpf(20) / 19) / 52,
        'complex_standardized_third_derivative_factor':
            48 * (mp.mpf(12) / 11)**4 / (mp.mpf(19) / 20)**mp.mpf('1.5'),
        'real_third_derivative_factor': 48 * (mp.mpf(256) / 255)**4,
        'real_fourth_derivative_factor': 768 * (mp.mpf(256) / 255)**5,
        'tilted_variance_deviation_bound': mp.mpf(50) / 1024,
        'normal_integral_lower_coefficient': (1 - mp.exp(-1)) / mp.sqrt(2 * mp.pi),
        'legendre_fourth_derivative_factor':
            800 * (mp.mpf(20) / 19)**4 + 7500 * (mp.mpf(20) / 19)**5,
        'maximum_log_corollary_relative_error': mp.mpf(5000) / 16384,
        'uniform_variance_error_factor': 4 + 2 * mp.exp(-2),
    }
    for key, bound in (
        ('tilted_esseen_constant', 200), ('tilt_path_radius_times_delta', mp.mpf(1)/48),
        ('complex_standardized_third_derivative_factor', 78), ('real_third_derivative_factor', 50),
        ('real_fourth_derivative_factor', 800), ('tilted_variance_deviation_bound', mp.mpf(1)/20),
        ('legendre_fourth_derivative_factor', 11000),
        ('maximum_log_corollary_relative_error', mp.mpf(1)/2),
    ):
        checks.upper('proof_constant_' + key, parameters[key], bound)
    checks.upper('proof_constant_normal_integral_lower_bound', mp.mpf(1)/4,
                 parameters['normal_integral_lower_coefficient'])
    checks.upper('proof_constant_cubic_remainder', mp.mpf(11000)/24, 500)
    for k in range(3, 31):
        checks.upper('classical_cumulant_condition_conversion', 2*k*(k-1), mp.mpf(12)**(k-2), k)
    gaussian_integral = mp.quad(lambda s: s*s*mp.exp(-s*s/4), [-mp.inf, 0, mp.inf])
    checks.identity('gaussian_smoothing_integral', gaussian_integral, 4 * mp.sqrt(mp.pi))
    return {key: decimal(value) for key, value in parameters.items()}


def binomial_logtail_beta(m, k, p):
    """Compute P(Bin(m,p)>=k) by a beta integral scaled at its endpoint.

    The logarithm of the beta integrand is concave. After u=p-y/L', its
    derivative at y=0 is -1, so the transformed integrand is <=exp(-y).
    Truncating at BETA_CUTOFF therefore omits at most exp(-BETA_CUTOFF)
    in the transformed integral. mpmath quadrature itself is not certified.
    """
    p = scalar(p)
    q = 1 - p
    assert 1 <= k <= m and k > m*p + q
    a, b = mp.mpf(k - 1), mp.mpf(m - k)
    derivative = (mp.mpf(k) - m*p - q) / (p*q)
    endpoint = p * derivative
    cutoff = min(mp.mpf(BETA_CUTOFF), endpoint)

    def integrand(y):
        if y >= endpoint:
            return mp.mpf(0)
        return mp.exp(a * mp.log1p(-y / endpoint)
                      + b * mp.log1p(y / (derivative*q)))

    intervals = sorted({mp.mpf(0), cutoff} | {mp.mpf(x) for x in (1, 4, 16, 64)
                                             if x < cutoff})
    integral, quadrature_error = mp.quad(integrand, intervals, error=True)
    assert integral > 0
    omitted_bound = mp.exp(-cutoff) if endpoint > cutoff else mp.mpf(0)
    log_prefactor = (mp.loggamma(m + 1) - mp.loggamma(k) - mp.loggamma(m - k + 1)
                     + a*mp.log(p) + b*mp.log(q) - mp.log(derivative))
    return log_prefactor + mp.log(integral), {
        'endpoint_derivative': decimal(derivative), 'transformed_upper_endpoint': decimal(endpoint),
        'quadrature_cutoff': decimal(cutoff), 'transformed_integral': decimal(integral),
        'quadrature_error_estimate_not_certified': decimal(quadrature_error),
        'analytic_omitted_integral_upper_bound': decimal(omitted_bound),
        'analytic_omitted_integral_bound_over_computed_integral': decimal(omitted_bound / integral),
    }


def exact_binomial_tail(m, k, p):
    p = Fraction(p)
    q = 1 - p
    return sum((comb(m, j) * p**j * q**(m-j) for j in range(k, m+1)), Fraction(0))


def log_normal_tail(x):
    return mp.log(mp.erfc(x / mp.sqrt(2)) / 2)


def verify_binomial_tails(checks):
    crosschecks, boundaries = [], []
    for p_string in ('1/3', '1/2', '2/5'):
        p = scalar(p_string)
        for m in (30, 90, 300):
            k = int(mp.ceil(m*p + 2*mp.sqrt(m*p*(1-p))))
            logtail, diagnostics = binomial_logtail_beta(m, k, p_string)
            exact = scalar(exact_binomial_tail(m, k, p_string))
            checks.identity('scaled_beta_vs_exact_binomial_pmf', logtail, mp.log(exact), (m, k, p_string))
            scipy_value = float(binom.logsf(k - 1, m, float(p)))
            double_error = abs(logtail - mp.mpf(scipy_value))
            assert double_error < mp.mpf(DOUBLE_TOLERANCE), (m, k, p_string, double_error)
            checks.counts['scaled_beta_vs_scipy_logsf'] += 1
            checks.maxima['scaled_beta_vs_scipy_logsf_absolute_error'] = max(
                checks.maxima['scaled_beta_vs_scipy_logsf_absolute_error'], double_error)
            crosschecks.append({'m': m, 'k': k, 'p': p_string, 'log_tail': decimal(logtail),
                                'scipy_absolute_error': decimal(double_error), **diagnostics})

    p = mp.mpf(1) / 3
    q = 1 - p
    for label, exponent, sizes in (
        ('gaussian_boundary', mp.mpf(1)/3, (10**6, 10**9, 10**12)),
        ('cubic_boundary', mp.mpf(1)/2, (10**8, 10**12, 10**16)),
    ):
        for m in sizes:
            delta = 1 / mp.sqrt(m*p*q)
            x = delta**(-exponent)
            threshold = m*p + x / delta
            k = int(mp.ceil(threshold))
            logtail, diagnostics = binomial_logtail_beta(m, k, '1/3')
            lognormal = log_normal_tail(x)
            r = p + p*q*delta*x
            assert 0 < r < 1
            I = m * (r*mp.log(r/p) + (1-r)*mp.log((1-r)/q))
            correction = x*x/2 - I
            h = mp.log(r*q / (p*(1-r))) / delta
            v = r*(1-r)/(p*q)
            a = h*mp.sqrt(v)
            log_J = a*a/2 + log_normal_tail(a)
            full_error = abs(mp.expm1(logtail - lognormal - correction))
            saddle_error = abs(mp.expm1(logtail + I - log_J))
            theorem_domain = x <= 1 / (2048*delta)
            log_corollary_domain = (delta <= mp.mpf(1)/16384
                                    and x <= 1/(16384*delta))
            if theorem_domain:
                checks.upper('binomial_full_correction_bound_in_theorem_domain',
                             full_error, 2500*delta*(1+x), (label, m))
                checks.upper('binomial_saddlepoint_bound_in_theorem_domain',
                             saddle_error, 1800*delta*(1+x), (label, m))
            kappa3 = (q-p)*delta
            corrected_log_ratio = logtail - lognormal - kappa3*x**3/6
            if label == 'gaussian_boundary':
                illustrated_ratio = mp.exp(logtail - lognormal)
                limit = mp.exp(mp.mpf(1)/18)
            else:
                illustrated_ratio = mp.exp(corrected_log_ratio)
                limit = mp.exp(-mp.mpf(1)/36)
            boundaries.append({
                'illustration': label, 'm': m, 'p': '1/3', 's': '1',
                'delta': decimal(delta), 'x': decimal(x), 'integer_count_threshold_ge': k,
                'real_count_threshold': decimal(threshold),
                'inside_full_correction_theorem_domain': bool(theorem_domain),
                'inside_log_corollary_domain': bool(log_corollary_domain),
                'interpretation': ('Finite theorem-domain check and asymptotic illustration.' if theorem_domain
                                   else 'Asymptotic illustration outside the conservative finite theorem domain.'),
                'log_tail': decimal(logtail), 'log_normal_tail': decimal(lognormal),
                'legendre_correction': decimal(correction),
                'full_correction_relative_error': decimal(full_error),
                'full_correction_error_bound': decimal(2500*delta*(1+x)),
                'saddlepoint_relative_error': decimal(saddle_error),
                'saddlepoint_error_bound': decimal(1800*delta*(1+x)),
                'cubic_corrected_log_ratio': decimal(corrected_log_ratio),
                'illustrated_ratio': decimal(illustrated_ratio),
                'proved_asymptotic_limit': decimal(limit),
                'absolute_difference_from_limit': decimal(abs(illustrated_ratio-limit)),
                **diagnostics,
            })
    return {'small_binomial_crosschecks': crosschecks, 'boundary_illustrations': boundaries}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path,
                        default=Path(__file__).resolve().parents[1] / 'reports' / 'quantitative_extensions_verification.json',
                        help='JSON report path (default: repository reports/quantitative_extensions_verification.json)')
    args = parser.parse_args()
    mp.mp.dps = PRECISION
    checks = Checks()
    catalog = [prepare_model(specification) for specification in model_catalog()]
    constants = verify_constants(checks)
    tilts = verify_tilts(catalog, checks)
    fixed = verify_fixed_samples(catalog, checks)
    thresholds = verify_envelopes_and_loads(catalog, checks)
    binomial = verify_binomial_tails(checks)
    script = Path(__file__).resolve()
    helper = script.with_name('verify_optimal_rate.py')
    report = {
        'description': 'Deterministic supplementary numerical and exact finite-state checks. '
                       'The proofs establish the theorems; these checks are not formal or interval-certified verification.',
        'status': 'passed',
        'metadata': {
            'python_version': platform.python_version(), 'python_implementation': platform.python_implementation(),
            'platform': platform.platform(), 'mpmath_version': version('mpmath'), 'scipy_version': version('scipy'),
            'mpmath_decimal_precision': PRECISION, 'high_precision_relative_or_absolute_tolerance': TOLERANCE,
            'scipy_logtail_absolute_tolerance': DOUBLE_TOLERANCE,
            'script_sha256': hashlib.sha256(script.read_bytes()).hexdigest(),
            'optimal_rate_helper_sha256': hashlib.sha256(helper.read_bytes()).hexdigest(),
            'manuscript_source_path': 'papers/quantitative_extensions.tex',
            'model_parameter_seed': None, 'model_selection': 'Explicit rational parameters; no random generator.',
            'fixed_sample_arithmetic': 'Exact fractions for the union-state probabilities, weights, mean, and variance.',
            'poisson_arithmetic': 'Exhaustive independent batch-presence states at 80 decimal digits.',
            'large_binomial_method': 'Endpoint-scaled incomplete-beta integral, truncated at y=256 when necessary.',
            'large_binomial_error_scope': 'The omitted-integral bound is analytic; mpmath quadrature and arithmetic are not certified.',
        },
        'case_totals': {
            'explicit_finite_models': len(catalog), 'poisson_model_time_cases': len(original_models()) * len(TIMES),
            'poisson_tilt_cases': len(tilts), 'fixed_sample_model_size_cases': len(fixed),
            'envelope_thresholds': len(thresholds['envelope_integrals']),
            'two_rate_envelope_cases': len(thresholds['two_rate_integrals']),
            'weighted_coordinate_load_cases': len(thresholds['weighted_loads']),
            'small_binomial_crosschecks': len(binomial['small_binomial_crosschecks']),
            'large_binomial_illustrations': len(binomial['boundary_illustrations']),
            'large_binomial_cases_inside_full_correction_domain': sum(
                int(row['inside_full_correction_theorem_domain'])
                for row in binomial['boundary_illustrations']),
        },
        'parameters': {'poisson_times': list(TIMES), 'fixed_sample_sizes': list(SAMPLE_SIZES),
                       'tilt_profiles_h_times_delta': list(TILT_PROFILES),
                       'thresholds_for_envelope_checks': list(ENVELOPE_THRESHOLDS),
                       'scaled_beta_cutoff': BETA_CUTOFF},
        'models': [{'name': m['name'], 'batches': m['batches'], 'raw_rate_weights': m['raw_rates'],
                    'species_weights': m['weight_strings'], 'positive_probability_batch_weight_bound': str(m['W'])}
                   for m in catalog],
        **checks.summary(),
        'proof_constants': constants, 'poisson_tilt_checks': tilts,
        'fixed_sample_checks': fixed, 'threshold_checks': thresholds,
        'binomial_tail_checks': binomial,
    }
    destination = args.output.expanduser().resolve()
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(dict(sorted(checks.counts.items())), indent=2))
    print('All quantitative-extension checks passed. Report:', destination)


if __name__ == '__main__':
    main()
