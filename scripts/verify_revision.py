#!/usr/bin/env python3
"""Deterministic checks of the revised inequalities (not proofs).

Requires NumPy, SciPy, SymPy and mpmath. Finite discovery laws are enumerated
exhaustively over batch-presence states, with fixed model parameters; no
Monte Carlo sampling is used. Floating-point checks supplement, rather than
establish, the manuscript's proofs.
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
import mpmath as mp
import sympy as sp
from scipy.special import logsumexp, gammaln


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
        'model_parameter_seed': None,
        'model_selection': 'Fixed, explicitly specified parameters; no random generator.',
        'mpmath_decimal_precision': mp.mp.dps,
    }


def enumerate_model(batches, probabilities, d, t):
    probs = np.asarray(probabilities, dtype=float)
    if len(batches) != len(probs) or np.any(probs <= 0) or not np.isclose(probs.sum(), 1):
        raise ValueError('The probabilities must be positive and sum to one.')
    inc = np.zeros((len(batches), d))
    for j, batch in enumerate(batches):
        if len(set(batch)) != len(batch) or any(i < 0 or i >= d for i in batch):
            raise ValueError('Invalid batch.')
        inc[j, list(batch)] = 1
    lam = probs @ inc
    if np.any(lam == 0):
        raise ValueError('Every species must have positive incidence.')
    pairs = inc.T @ (probs[:, None] * inc)
    kappa = pairs.sum(axis=1)/lam
    sizes = inc.sum(axis=1).astype(int)
    B = int(sizes.max())
    Bi = np.max(sizes[:, None] * inc, axis=0)
    q = -np.expm1(-t*probs)
    states = np.array(list(itertools.product((0, 1), repeat=len(batches))))
    weights = np.prod(np.where(states, q, 1-q), axis=1)
    seen = (states @ inc > 0).astype(float)
    assert abs(weights.sum()-1) < 1e-12
    score = seen.sum(axis=1)
    mean = weights @ score
    V = weights @ (score-mean)**2
    vi = np.exp(-t*lam)*(-np.expm1(-t*lam))
    return weights, seen, score, mean, V, vi, kappa, B, Bi


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path,
                        default=Path(__file__).resolve().parents[1]/'reports'/'revision_verification_results.json',
                        help='JSON report path (default: repository reports/revision_verification_results.json)')
    args = parser.parse_args()
    models = [
        ('singletons', [(0,), (1,), (2,), (3,)], [.05,.15,.3,.5], 4),
        ('pairs', [(0,1),(2,3)], [.25,.75], 4),
        ('star', [(0,1),(0,2),(0,3)], [.1,.3,.6], 4),
        ('mixed', [(),(0,),(1,),(0,1),(1,2,3),(0,2)], [.05,.1,.2,.15,.3,.2], 4),
        ('heterogeneous', [(0,), (1,2), (3,4,5)], [.5,.499,.001], 6),
        ('overlapping', [(0,1,2),(1,3),(0,4),(2,4),(4,5),(5,)], [.1,.2,.15,.1,.25,.2], 6),
    ]
    totals = {'finite_models':len(models), 'signed_variance_checks':0,
              'fractional_mgf_checks':0, 'actual_variance_mgf_checks':0,
              'coefficient_load_checks':0, 'threshold_envelope_checks':0}
    records=[]
    max_violation = 0.0
    for name,batches,probs,d in models:
        for t in (.02,.2,1.,4.,12.,50.):
            weights,seen,score,m,V,vi,kappa,B,Bi = enumerate_model(batches,probs,d,t)
            logweights = np.full(weights.shape, -np.inf)
            np.log(weights, out=logweights, where=weights>0)
            vectors = [np.ones(d), np.array([(-1)**i*(i+1) for i in range(d)]),
                       np.arange(d)/max(1,d-1)]
            for w in vectors:
                z = seen @ w
                true = weights @ (z-weights @ z)**2
                bound = np.sum(kappa*w*w*vi)
                assert true <= bound + 2e-11, (name,t,true,bound)
                max_violation = max(max_violation, float(true-bound))
                totals['signed_variance_checks'] += 1
            for theta in (-1.,-.2,.2,1.):
                logmgf = logsumexp(logweights + theta*(score-m))
                fractional = np.sum(vi/Bi * (np.expm1(Bi*abs(theta))-Bi*abs(theta)))
                global_bound = vi.sum()/B * (np.expm1(B*abs(theta))-B*abs(theta))
                assert logmgf <= fractional + 2e-11, (name,t,theta)
                assert fractional <= global_bound + 2e-11
                totals['fractional_mgf_checks'] += 1
            for a in (.02,.1,.24,.45):
                for sign in (-1,1):
                    theta=sign*a/B
                    logmgf=logsumexp(logweights+theta*(score-m))
                    bound=V*theta*theta/(2*(1-2*a))
                    assert logmgf <= bound+2e-11, (name,t,theta,logmgf,bound)
                    totals['actual_variance_mgf_checks']+=1
            records.append({'model':name,'t':t,'variance':float(V),
                            'species_variance_bound':float(np.sum(kappa*vi)),
                            'global_variance_bound':float(B*vi.sum())})
    a=sp.symbols('a')
    q=1+a
    difference=sp.expand(1-q*a/3-q*q*a/2-q*(1-2*a))
    assert sp.simplify(difference-(a/6+2*a*a/3-a**3/2))==0
    holder_grid = np.linspace(.0001,.4999,5000)
    for x in holder_grid:
        qx=1+x
        D=1-qx*x/3-qx*qx*x/2
        assert D>0 and qx/D <= 1/(1-2*x) + 1e-10
        z=qx*qx*x/(1-qx*x/3)
        assert z<2
    totals['holder_parameter_checks']=len(holder_grid)
    loads=[]
    for B in range(2,8):
        L=B*(4*B-3)/2
        for n in (1,3,20,1000,100000):
            load=B/2+2*B*n*(B-1)/(n+2)+B*(B-1)/(2*(n+1))
            assert load<=L+1e-12
            loads.append({'B':B,'N':n,'load':load,'L_B':L})
            totals['coefficient_load_checks']+=1
    for r in (2,3,4,5,10,20,50,100):
        for k in range(r):
            x=np.unique(np.r_[0.,np.linspace(.001,3*r+5,301),float(k)])
            positive=x>0
            logpmf=-x[positive]+k*np.log(x[positive])-gammaln(k+1)
            # The better of the old envelope and the new linear-r envelope
            assert np.all(logpmf <= 1-x[positive]/r + 2e-12)
            assert np.all(logpmf <= (r-1)*np.log(2)-x[positive]/2 + 2e-12)
            totals['threshold_envelope_checks']+=2*len(logpmf)
    mp.mp.dps=70
    p=-mp.expm1(-1); q=p*(1-p)
    blocks=[]
    for K in range(2,11):
        var_ratio=mp.mpf(0); prob_tail=mp.mpf(0); mean_tail=mp.mpf(0)
        for k in range(1,K+100):
            delta=k*k-K*K
            rate=mp.power(2,-2*delta-(k-K))
            discovery=-mp.expm1(-rate)
            var_ratio+=mp.power(2,2*delta)*mp.exp(-rate)*discovery
            if k>K:
                prob_tail+=rate
                mean_tail+=mp.power(2,delta)*discovery
        blocks.append({'K':K,'variance_ratio':float(var_ratio),
                       'limit_variance_ratio':float(1+q),
                       'later_occurrence_upper_bound':float(prob_tail),
                       'later_mean_over_block_size':float(mean_tail)})
    assert abs(blocks[-1]['variance_ratio']-float(1+q))<1e-10
    assert blocks[-1]['later_occurrence_upper_bound']<1e-10
    totals['counterexample_parameter_checks']=len(blocks)
    out={'description':
         'Deterministic supplementary checks with fixed model parameters: exhaustive '
         'batch-presence enumeration, symbolic algebra, and finite numerical grids; '
         'no Monte Carlo sampling. These checks are not proof or formal verification.',
         'metadata':run_metadata(),
         'totals':totals, 'maximum_signed_variance_violation':max_violation,
         'holder_polynomial_identity':str(difference), 'finite_models':records,
         'load_sharpness':loads, 'unbounded_counterexample':blocks}
    dest=args.output.expanduser().resolve()
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(totals,indent=2)); print('All revision checks passed. Report:',dest)

if __name__=='__main__':
    main()
