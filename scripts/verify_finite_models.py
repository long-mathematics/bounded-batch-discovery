#!/usr/bin/env python3
"""Deterministic finite-state checks for the bounded-batch discovery draft.

Requires Python 3, NumPy, and SciPy. These checks are not mathematical proofs.
Two model parameter vectors are selected deterministically with a fixed
pseudo-random seed. No Monte Carlo sampling is used: probabilities and bracket
moments are solved from the finite-state Kolmogorov/Feynman--Kac equations.
"""
from __future__ import annotations
import argparse
import hashlib
from importlib.metadata import PackageNotFoundError, version
import itertools
import json
from pathlib import Path
import platform
import sys

if not __debug__:
    raise SystemExit('Verification requires assertions; run Python without -O or -OO.')

import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import gammaln, logsumexp
from scipy.stats import poisson


MODEL_SEED = 20261112
THRESHOLDS = (1, 2, 3)
OBSERVATION_TIMES = (.2, 1.0, 4.0, 12.0)


def run_metadata():
    dependencies = {}
    for package in ('numpy', 'scipy', 'sympy', 'mpmath'):
        try:
            dependencies[package] = version(package)
        except PackageNotFoundError:
            dependencies[package] = None
    return {
        'python_version': platform.python_version(),
        'python_implementation': platform.python_implementation(),
        'python_build': sys.version,
        'platform': platform.platform(),
        'package_versions': dependencies,
        'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'model_parameter_seed': MODEL_SEED,
        'random_generator': 'numpy.random.default_rng (PCG64)',
        'thresholds': list(THRESHOLDS),
        'observation_times': list(OBSERVATION_TIMES),
        'ode_solver': {'method': 'RK45', 'rtol': 2e-10, 'atol': 2e-12},
    }


def checked_model(batches: list[tuple[int, ...]], probs: np.ndarray, d: int):
    probs = np.asarray(probs, dtype=float)
    if len(batches) != len(probs) or np.any(probs < 0):
        raise ValueError('Invalid batch probabilities')
    if not np.isclose(probs.sum(), 1.0):
        raise ValueError('Probabilities must sum to one')
    inc = np.zeros((len(batches), d))
    for j, batch in enumerate(batches):
        if len(set(batch)) != len(batch) or any(i < 0 or i >= d for i in batch):
            raise ValueError('Invalid species set')
        inc[j, list(batch)] = 1.0
    lam = probs @ inc
    pair = inc.T @ (probs[:, None] * inc)
    B = max(map(len, batches))
    return probs, inc, lam, pair, B


def state_audit(batches, probs, d, r, t):
    probs, inc, lam, pair, B = checked_model(batches, probs, d)
    states = np.array(list(itertools.product(range(r+1), repeat=d)), dtype=int)
    lookup = {tuple(s): k for k, s in enumerate(states)}
    count = len(states)
    generator = np.zeros((count, count))
    for k, state in enumerate(states):
        for batch, prob in zip(batches, probs):
            nxt = state.copy()
            nxt[list(batch)] = np.minimum(r, nxt[list(batch)] + 1)
            j = lookup[tuple(nxt)]
            generator[j, k] += prob
        generator[k, k] -= 1.0
    score = (states >= r).sum(axis=1).astype(float)

    def rhs(s, y):
        p, u, w = y.reshape(3, count)
        x = lam[None, :] * max(t-s, 0.0)
        need = r-1-states
        h = np.zeros_like(x + states, dtype=float)
        for i in range(d):
            eligible = need[:, i] >= 0
            k = need[eligible, i]
            z = float(x[0, i])
            if z == 0:
                h[eligible, i] = (k == 0)
            else:
                h[eligible, i] = np.exp(-z + k*np.log(z)-gammaln(k+1))
        jumps = h @ inc.T
        q = (jumps*jumps) @ probs
        return np.concatenate((generator @ p,
                               generator @ u + q*p,
                               generator @ w + 2*q*u))

    init = np.zeros(3*count)
    init[lookup[tuple([0]*d)]] = 1.0
    sol = solve_ivp(rhs, (0, t), init, rtol=2e-10, atol=2e-12)
    if not sol.success:
        raise RuntimeError(sol.message)
    p, u, w = sol.y[:, -1].reshape(3, count)
    mean = p @ score
    V = max(float(p @ (score*score)-mean*mean), 0.0)
    EQ = float(u.sum())
    varQ = max(float(w.sum()-EQ*EQ), 0.0)
    envelope = 1.0 if r == 1 else min(2**(2*r-1), np.exp(2)*r)
    L = B*(4*B-3)/2 * envelope
    assert abs(EQ-V) < 3e-7*(1+V), (r,t,EQ,V)
    assert varQ <= L*V + 3e-7, (r,t,varQ,L*V)
    marginal_probabilities = poisson.sf(r-1, t*lam)
    exact_mean = float(marginal_probabilities.sum())
    marginal_variance = float(np.sum(marginal_probabilities*(1-marginal_probabilities)))
    mean_error = abs(float(mean)-exact_mean)
    assert mean_error < 3e-7*(1+exact_mean), (r,t,mean,exact_mean)
    assert V >= marginal_variance-3e-7, (r,t,V,marginal_variance)
    assert V <= B*marginal_variance+3e-7, (r,t,V,B*marginal_variance)
    mgf_checks = 0
    if r == 1:
        x = np.exp(-t*lam)
        v0 = float(np.sum(x*(1-x)))
        exactV = float(np.sum(np.exp(-t*(lam[:,None]+lam[None,:]-pair))
                              -np.exp(-t*(lam[:,None]+lam[None,:]))))
        assert abs(exactV-V) < 3e-7*(1+V)
        assert V >= v0-3e-7 and V <= B*v0+3e-7
        for theta in (-1.0, -0.25, 0.25, 1.0):
            logmgf = logsumexp(np.log(np.maximum(p, 1e-300)) + theta*(score-mean))
            upper = v0/B * (np.expm1(B*abs(theta))-B*abs(theta))
            assert logmgf <= upper + 3e-7, (theta,logmgf,upper)
            mgf_checks += 1
    return {'r':r,'time':t,'mean':float(mean),'variance':V,
            'mean_bracket':EQ,'variance_bracket':varQ,
            'bracket_bound':float(L*V),'bracket_identity_error':abs(EQ-V),
            'threshold_envelope_constant':float(envelope),
            'bracket_load_constant':float(L),
            'exact_threshold_mean':exact_mean,
            'threshold_mean_identity_error':mean_error,
            'marginal_variance_sum':marginal_variance,
            'variance_sandwich_lower_bound':marginal_variance,
            'variance_sandwich_upper_bound':float(B*marginal_variance),
            'bernoulli_mgf_checks':mgf_checks,
            'state_count':count}


def fixed_sample_audit(batches, probs, d, n):
    probs, inc, lam, c, B = checked_model(batches, probs, d)
    masks = [sum(1 << i for i in a) for a in batches]
    p = np.zeros(1 << d); p[0] = 1
    for _ in range(n):
        nxt = np.zeros_like(p)
        for s in range(1 << d):
            for mask, pr in zip(masks, probs):
                nxt[s | mask] += p[s]*pr
        p = nxt
    score = np.array([s.bit_count() for s in range(1 << d)])
    mean = p @ score
    Vfixed = float(p @ (score*score)-mean**2)
    a = lam[:, None]; b = lam[None, :]
    u = a+b-c; w = a+b-a*b
    fixed_formula = float(np.sum(np.maximum(1-u,0)**n
                        -(1-a)**n*(1-b)**n))
    Vpois = float(np.sum(np.exp(-n*u)-np.exp(-n*(a+b))))
    delta = lambda x: np.exp(-n*x)-np.maximum(1-x,0)**n
    E = float(np.sum(delta(w)-delta(u)))
    R = float(np.sum(np.exp(-n*w)-np.exp(-n*(a+b))))
    f = float(np.sum(lam*np.exp(-n*lam)))
    rem_bound = 0.5*(n*np.sum(lam**2*np.exp(-n*lam/2)))**2
    assert abs(Vfixed-fixed_formula) < 1e-8
    assert abs(Vfixed-(Vpois+E-R)) < 1e-8
    assert -1e-10 <= R-n*f*f <= rem_bound+1e-10
    return {'n':n,'variance_fixed':Vfixed,'variance_poisson':Vpois,
            'correction':n*f*f,'absolute_remainder':Vfixed-Vpois+n*f*f,
            'kernel_identity_error':abs(Vfixed-(Vpois+E-R))}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path,
                        default=Path(__file__).resolve().parents[1]/'reports'/'verification_results.json',
                        help='JSON report path (default: repository reports/verification_results.json)')
    args = parser.parse_args()
    rng = np.random.default_rng(MODEL_SEED)
    models = [
        ('singletons',4,[(i,) for i in range(4)],np.array([.05,.15,.30,.50])),
        ('disjoint_pairs',4,[(0,1),(2,3)],np.array([.25,.75])),
        ('star',4,[(0,1),(0,2),(0,3)],np.array([.1,.3,.6])),
        ('overlapping_pairs',4,list(itertools.combinations(range(4),2)),
         rng.dirichlet(np.ones(6))),
        ('mixed_batches',4,[(),(0,),(1,),(0,1),(1,2,3),(0,2)],
         rng.dirichlet(np.ones(6))),
    ]
    results = {'description':
               'Deterministic finite-state numerical sanity checks, not proofs: '
               'fixed-seed pseudo-random model parameters, Kolmogorov/Feynman--Kac '
               'ODE integration, and exact finite-sample state enumeration; no Monte Carlo sampling.',
               'metadata':run_metadata(),
               'models':{}}
    for name,d,batches,probs in models:
        checks=[]
        for r in THRESHOLDS:
            for t in OBSERVATION_TIMES:
                checks.append(state_audit(batches,probs,d,r,t))
        fixed=[fixed_sample_audit(batches,probs,d,n) for n in (1,2,5,20,80)]
        results['models'][name]={'species_count':d, 'batches':batches,
                                'probabilities':probs.tolist(),
                                'bracket_checks':checks,'fixed_sample_checks':fixed}
    state_checks = [z for model in results['models'].values() for z in model['bracket_checks']]
    results['totals']={
        'finite_models':len(results['models']),
        'bracket_checks':len(state_checks),
        'fixed_sample_checks':sum(len(model['fixed_sample_checks']) for model in results['models'].values()),
        'bernoulli_mgf_checks':sum(z['bernoulli_mgf_checks'] for z in state_checks),
        'threshold_mean_checks':len(state_checks),
        'threshold_variance_sandwich_checks':len(state_checks),
        'maximum_bracket_identity_error':max(z['bracket_identity_error'] for z in state_checks),
        'maximum_threshold_mean_identity_error':max(z['threshold_mean_identity_error'] for z in state_checks)}
    dest=args.output.expanduser().resolve()
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(json.dumps(results,indent=2)+'\n')
    print(json.dumps(results['totals'],indent=2))
    print('All numerical checks passed. Report:',dest)

if __name__ == '__main__':
    main()
