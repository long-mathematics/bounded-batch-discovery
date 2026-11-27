# Research status and remaining questions

**Manuscript revision: November 26, 2026.**

The uniform Poissonized first-discovery Berry–Esseen target is proved in the
optimal-rate manuscript. The quantitative-extensions manuscript adds relative
tail approximations, uniform fixed-sample moment corrections, and sharper
higher-coverage estimates. The fixed-sample manuscript now proves a zero-free
MGF disk and the rate $CW\log(n+2)/s_{h,n}$, with a signed-weight extension.
The optimal fixed-sample rate is proved in the regime
$s_{h,n}^2\ge\eta nW^2$ with a constant depending on $\eta>0$.
The full uniform optimal fixed-sample and higher-coverage Berry–Esseen
bounds remain unresolved in this repository.

Here, **proved** means that a proof is supplied in the cited research draft.
The internal proof audits and deterministic calculations do not constitute
independent peer review or formal verification. **Open here** describes the
scope of these manuscripts; it does not assert the absence of results for
special cases elsewhere in the literature.

## Manuscripts and notation

| Reference | Manuscript | Principal role |
|---|---|---|
| **F** | [A fluctuation theory for bounded-batch species discovery](../papers/bounded_batch_discovery.tex) · [PDF](../papers/bounded_batch_discovery.pdf) | Fixed-law CLT, variance comparisons, concentration, bracket estimates, de-Poissonization, and extensions. |
| **O** | [Optimal Berry–Esseen bounds for bounded-batch discovery](../papers/optimal_berry_esseen.tex) · [PDF](../papers/optimal_berry_esseen.pdf) | Zero-free MGF, variance-linear cumulants, and the optimal-order Poissonized first-discovery bound. |
| **Q** | [Further quantitative bounds for bounded-batch discovery](../papers/quantitative_extensions.tex) · [PDF](../papers/quantitative_extensions.pdf) | Relative tails and their sharp windows; uniform fixed-sample mean and variance corrections; improved threshold envelopes. |
| **S** | [Fixed-sample zero-free bounds for bounded-batch discovery](../papers/fixed_sample_bounds.tex) · [PDF](../papers/fixed_sample_bounds.pdf) | Fixed-sample zero-free disk, rate with a logarithmic factor, Hoeffding-kernel operator bounds, and a sector counterexample. |

The original foundation draft is dated November 12, 2026, with its repository
revision dated November 26. The three companion manuscripts are dated November 26,
2026. Historical source bundles and audit records retain their own dates.

Let $D_n$ be the number of species discovered after $n$ independent batches,
and let $X(t)=D_{N(t)}$, where $N$ is an independent rate-one Poisson process.
Write $V(t)=\operatorname{Var}X(t)$ and $s_n^2=\operatorname{Var}D_n$.
For nonnegative weights $h_i$, the corresponding counts are $D_{h,n}$ and
$X_h(t)$, and the per-batch assumption is

$$
\sum_{i\in A_1}h_i\le W\quad\text{almost surely}.
$$

The unweighted bounded-cardinality model has $h_i=1$ and $W=B$. Bounded total
weight allows unbounded finite batch cardinalities. Unless explicitly stated
otherwise, the law is deterministic, the species space is countable, and
overlap is arbitrary.

## Status at a glance

| Question | Status | Scope and qualification |
|---|---|---|
| Poissonized first-discovery rate $CB/\sqrt V$ | **Proved in O.** | Constant 150; nonnegative weighted version $150W/\sqrt{V_h}$. Uniform in law, time, and overlap. The order is sharp. |
| Variance-linear first-discovery cumulants | **Proved in O.** | $\lvert\kappa_k(X_h)\rvert\le2k!V_h(4W)^{k-2}$ for $k\ge3$. |
| Gaussian relative tails and Cramér corrections | **Proved in Q.** | Poissonized first discovery; both signs and endpoint conventions. The general cumulant-to-tail transfer is classical. |
| Fixed-law, fixed-sample CLT | **Proved in F.** | Bounded batches and diverging actual fixed-sample variance; this is not a uniform changing-law Berry–Esseen theorem. |
| Uniform fixed-sample mean and corrected variance | **Proved in Q.** | Absolute errors $O(W)$ and $O(W^2)$ for every law and $n$; these error orders are sharp. |
| Fixed-sample zero-free MGF and rate $CW\log(n+2)/s_{h,n}$ | **Proved in S.** | Arbitrary overlap, independent nonidentical laws, countable species; signed weights with per-batch absolute total at most $W$. |
| Fixed-sample rate $CW/s_{h,n}$ under arbitrary overlap | **Open in full generality.** | Proved in S when $s_{h,n}^2\ge\eta nW^2$, with constant $C_\eta$. Removing the logarithm uniformly as the variance ratio decreases remains unresolved. |
| Fixed-sample Hoeffding-kernel flattening bounds | **Proved in S.** | Iid nonnegative weights; bounds are summable across orders with any fixed polynomial weight. Mixed-contraction control remains open. |
| Universal full zero-free PGF sector | **False, already for $B=2$.** | Two-community construction in S; the untilted variance tends to zero, so this does not disprove the optimal rate. |
| Higher-coverage bracket, CLT, and concentration bounds | **Improved in Q.** | The threshold parameter decreases from order $r$ to order $\sqrt r$. |
| Higher-coverage rate $CW/\sqrt{V_r}$, uniform in thresholds | **Open here.** | No proof or counterexample is supplied. An optimal-in-variance rate with a useful threshold factor also remains open in the general overlap model. |
| Signed weights and exact-frequency CLTs under variance divergence alone | **Open here.** | S supplies a fixed-sample signed-weight rate with a logarithmic factor. It does not settle the variance-only criterion or exact-frequency counts. |
| Unrestricted span-one local limits | **False in general.** | Disjoint pairs give even-valued counts. A suitable arithmetic version remains open here. |
| Joint occupancy-spectrum and functional CLTs | **Open here.** | Scalar limits, temporal covariances, and increment bounds are inputs, not joint or process convergence. |
| Finite second moment alone for unbounded batches | **Insufficient.** | F constructs a fixed-law non-Gaussian subsequential limit despite variance divergence. Additional Lindeberg or weight control is supplied for positive results. |

## 1. What the optimal first-discovery result establishes

Manuscript O proves, for every $V_h(t)>0$,

$$
d_K\!\left(\frac{X_h(t)-\mathbb EX_h(t)}{\sqrt{V_h(t)}},N(0,1)\right)
\le\min\left\{1,\frac{150W}{\sqrt{V_h(t)}}\right\}.
$$

The actual coefficient obtained in the proof is
$36/\sqrt\pi+864/(\pi\sqrt{2\pi})<131$; the theorem rounds it to 150.
Disjoint batches give $X=B\operatorname{Bin}(m,1/2)$, whose largest atom
forces a Kolmogorov error of order $B/\sqrt V$ against any continuous
reference distribution. The numerical constant is not claimed optimal.

The proof combines an exact bipartite independence-polynomial representation
with the actual-variance real MGF estimate from F. It permits unbounded
species dependency degree. It does not rely on the external dependency-graph
theorem discussed in [the earlier comparison note](berry_esseen.md).
The zero-free disk and the real MGF bound are both needed for the
variance-linear cumulant estimate; zero-freeness alone does not provide it.

## 2. Relative tails: a completed extension with classical provenance

For Poissonized nonnegative weighted first discovery, set

$$
S=\frac{X_h-\mathbb EX_h}{\sqrt V},\qquad
\delta=\frac W{\sqrt V},\qquad \overline\Phi=1-\Phi.
$$

For a sign $\epsilon\in\{-1,1\}$, put
$K_\epsilon(z)=\log\mathbb E e^{z\epsilon S}$.
For $0\le x\le(2048\delta)^{-1}$, Q constructs the unique small nonnegative
solution $h$ of $K_\epsilon'(h)=x$ and defines

$$
I_\epsilon(x)=hx-K_\epsilon(h),\qquad
\Lambda_\epsilon(x)=\frac{x^2}{2}-I_\epsilon(x).
$$

For either $\diamond\in\{>,\ge\}$, it proves the full-correction estimate

$$
\left|
\frac{\Pr(\epsilon S\diamond x)}
{\overline\Phi(x)e^{\Lambda_\epsilon(x)}}-1
\right|
\le2500\delta(1+x).
$$

The saddlepoint form is also supplied. In a stated smaller universal range,

$$
\log\frac{\Pr(\epsilon S\diamond x)}{\overline\Phi(x)}
=\frac{\epsilon\kappa_3(S)x^3}{6}
 +O\!\left(\delta(1+x)+\delta^2x^4\right).
$$

As $\delta\to0$, the resulting uniform windows are:

| Approximation with relative error tending to zero | Window |
|---|---|
| Gaussian tail $\overline\Phi(x)$ | $x=o(\delta^{-1/3})$ |
| Cubic-corrected tail $\overline\Phi(x)\exp\{\epsilon\kappa_3(S)x^3/6\}$ | $x=o(\delta^{-1/2})$ |
| Full correction $\overline\Phi(x)e^{\Lambda_\epsilon(x)}$ | $x=o(\delta^{-1})$ |

Poissonized singleton examples prove that the first two window exponents
cannot be enlarged uniformly over this model class. No optimality claim is
made for the numerical constants or the full-correction window. Endpoint
atoms are controlled explicitly; these tail results require neither a
density nor a specified lattice span. They do not establish a local limit
theorem or transfer relative tails to a deterministic sample size. See Q,
Sections 1–6.

**Attribution.** The transfer from cumulant bounds to corrected relative
tails is classical Cramér–Petrov theory. In particular,
[Döring–Jansen–Schubert, Theorem 2.1 and Section 6][DJS] already give the
general transfer needed for the order of these conclusions. The input from O
satisfies their cumulant condition with $\Delta=(48\delta)^{-1}$, since

$$
2k!(4\delta)^{k-2}\le(k-2)!(48\delta)^{k-2},\qquad k\ge3.
$$

Q applies that theory to bounded-batch discovery and supplies a
self-contained tilt proof of its conservative constants, endpoint handling,
cubic remainder, and model-specific window examples. The general transfer
principle is not presented as a new theorem of this project.

## 3. Fixed samples: moment corrections and a rate with a logarithmic factor

Let $b_{h,n}=\mathbb ED_{h,n}$, $m_h(t)=\mathbb EX_h(t)$, and

$$
f_h(t)=m_h'(t)=\sum_i h_i\lambda_i e^{-t\lambda_i},\qquad
\tau_{h,n}^2=V_h(n)-n f_h(n)^2.
$$

Q, Section 7, proves for every $n\ge1$ and every sampling law

$$
0\le b_{h,n}-m_h(n)\le2W,
\qquad
|s_{h,n}^2-\tau_{h,n}^2|\le(4+2e^{-2})W^2.
$$

The exact residual interpretation is

$$
\tau_{h,n}^2
=\operatorname{Var}\bigl(X_h(n)-f_h(n)N(n)\bigr)
=\min_{a\in\mathbb R}\operatorname{Var}\bigl(X_h(n)-aN(n)\bigr).
$$

It also has a positive decomposition into a residual first Poisson-chaos
variance and a nonnegative higher-order sum. Thus
$s_{h,n}^2\sim\tau_{h,n}^2$ whenever either variance diverges relative to
$W^2$, uniformly over changing laws. The sharper moment-dependent bounds in
Q require only a finite second moment of total batch weight; that extension
concerns moment corrections, not a CLT for arbitrary unbounded batches.

The uniform error orders are necessary. For uniform singleton sampling from
$n$ species with $n$ draws,

$$
b_n-m(n)\longrightarrow\frac1{2e},\qquad
s_n^2-\tau_n^2\longrightarrow\frac{3e^{-2}-e^{-1}}2>0.
$$

Disjoint-block or weight scaling multiplies these discrepancies by $B,B^2$
or $W,W^2$. This proves sharpness of the correction orders, not of a
fixed-sample normal-approximation rate. The classical singleton results of
[Hwang–Janson][HJ], including their uniform moment correction and local
limit theorem, are relevant benchmarks; the new correction covers arbitrary
overlap and nonnegative batch weights.

### A variance-cancellation example

Take $B=1$ and, for sample size $n$, use the uniform singleton law on
$M_n=\lceil n^{3/2}\rceil$ species. With $a_n=n/M_n$,

$$
V(n)=M_n e^{-a_n}(1-e^{-a_n})\sim n,
\qquad f(n)=e^{-a_n},
$$

whereas

$$
\tau_n^2=M_n e^{-2a_n}(e^{a_n}-1-a_n)
\sim\frac{\sqrt n}{2},
\qquad s_n^2\sim\frac{\sqrt n}{2}.
$$

Both variances diverge, but $V(n)/s_n^2\sim2\sqrt n\to\infty$.
This disproves uniform substitution of the Poissonized variance over
changing laws. It does not contradict the fixed-law variance equivalence
in F and does not disprove a fixed-sample bound $CB/s_n$.

### A distributional theorem now proved

S, Theorems 1.1–1.2, assumes only independent draws and
$\sum_{i\in A_k}|h_i|\le W$ and proves

$$
\mathbb Ee^{zD_{h,n}}\ne0\quad (|z|<1/(4W)),
\qquad
d_K\!\left(\frac{D_{h,n}-\mathbb ED_{h,n}}{s_{h,n}},Z\right)
\le\min\{1,CW\log(n+2)/s_{h,n}\}.
$$

The graph representation has an independent species side and a disjoint
union of cliques on the batch-choice side, one clique per draw. Its
occupation-ratio induction uses an arcsine angle budget of $\pi/4$ per
choice. Normalizing activities $M\mu_k(C)$ and letting $M\to\infty$
recovers the fixed-sample absence transform. Hurwitz's theorem handles
this limit and the countable-species exhaustion.

The rate is an application of [Michelen–Sahasrabudhe, Lemma 4.1 and
Corollary 9.3](https://doi.org/10.4171/JEMS/1530). A Harnack bound for
$\log\mathbb Ee^{zD/W}$ supplies the boundary control on their truncated
sector. This yields $\log(n+2)$, independently of cardinality, lattice
mesh, or the smallest weight. Signed weights follow by translation of
$D/W\in[-n,n]$ into $[0,2n]$.

### The remaining optimal-rate argument

The target is a uniform bound of order $CW/s_{h,n}$ for
$(D_{h,n}-b_{h,n})/s_{h,n}$, with no logarithm. S, Proposition 7.1,
proves that it would follow from the local actual-variance estimate

$$
\log\mathbb Ee^{\theta(D_{h,n}-\mathbb ED_{h,n})}
\le A s_{h,n}^2\theta^2,\qquad |\theta|\le a/W,
$$

with absolute $A,a>0$. The zero-free disk is now available, so such a
real-axis estimate would give variance-linear cumulants by a harmonic
coefficient bound and then the optimal Kolmogorov rate by smoothing.
The estimate has not been proved for the full class. Bounded differences
does prove it when $s_{h,n}^2\ge\eta nW^2$; S, Corollary 7.2, gives
$C_\eta W/s_{h,n}$ in that regime. Its dependence on $\eta$ cannot be
discarded.

For iid nonnegative weights, write $p_i=\Pr(i\in A)$,
$q_i=1-p_i$, and $\xi_i(A)=1_{\{i\in A\}}-p_i$. S, Proposition 5.1,
gives the canonical kernels and their operator bounds:

$$
H_k=(-1)^{k+1}\sum_i h_iq_i^{n-k}\xi_i^{\otimes k},
\qquad
\sqrt{\binom nk}\|H_k\|_{r\to k-r}\le W
\quad(1\le r<k).
$$

After excluding $p_i>1/2$ species, whose total weight is at most $2W$,
the right side improves to $W(2/3)^{k/2}$. The high-probability part is
at most $2W2^{-n/2}$ by a separate rank-one estimate. Thus the full
operator bounds are summable over $k$, also with a fixed polynomial
weight. These are operator norms, not the Hilbert–Schmidt norms in the
variance sum. A variance-sensitive estimate for mixed contractions
across orders is still missing.

The available scalar Poissonized theorem
normalizes $X_h(n)$ by $\sqrt{V_h(n)}$. It does not approximate the
residual $X_h(n)-f_h(n)N(n)$ at scale $\tau_{h,n}$, and does not provide
the joint local control needed to condition on $N(n)=n$. A joint
Kolmogorov approximation alone would not automatically control that
conditioning event, whose probability is of order $n^{-1/2}$.

The residual jump is $H_C-f_h(n)$ and can have either sign. Its squared
bracket introduces signed cross terms, so the positive bracket
decomposition used for the actual-variance MGF estimate does not simply
transfer. The fixed-sample absence probabilities
$[1-\Pr(A_1\cap S\ne\varnothing)]^n$ lack the independent-clock
factorization used in O, but the new clique representation now handles
their zero-freeness. An adequate residual Stein estimate, a joint
conditional Fourier argument, or a direct actual-variance cumulant bound
could remove the remaining logarithm. Moment corrections and a qualitative coupling do not settle
fixed-sample relative tails or moderate deviations either.

### A proved obstruction to the global sector route

For two communities of $M\ge3$ species, choose a community with
probability $1/2$ and then a uniform pair inside it. With $B=2$, S proves

$$
2^{2n}\mathbb E(2^{-n/M}z)^{D_n}\longrightarrow z^{2M}+2z^M.
$$

Thus there are PGF roots asymptotic to
$2^{-(n-1)/M}e^{i\pi/M}$. Their angles can be arbitrarily small while
$B$ stays two. No universal global sector $|\arg w|<c/B$ works, so the
sector theorem of Michelen–Sahasrabudhe cannot remove the logarithm by
that route. Here $s_n^2\sim2M(1-1/M)^n\to0$; this is not a
counterexample to the proposed optimal Berry–Esseen inequality.

## 4. Higher coverage: a better threshold parameter and an open optimal rate

For integer thresholds $1\le r_i\le r$, let

$$
Y_r(t)=\sum_i h_i1_{\{N_i(t)\ge r_i\}},\qquad
V_r(t)=\operatorname{Var}Y_r(t).
$$

Q, Section 8, proves the exact integrated mode envelope

$$
\ell_r=1+\sum_{k=1}^{r-1}\frac{e^{-k}k^k}{k!}
=\int_0^\infty\max_{0\le k<r}\Pr(\operatorname{Poi}(x)=k)\,dx
\sim\sqrt{\frac{2r}{\pi}}.
$$

It supplies a positive bracket decomposition with load at most
$2\ell_rW^2$. With $\eta_r=\ell_rW^2/V_r$, the generic martingale bounds
therefore give

$$
d_K\!\left(\frac{Y_r-\mathbb EY_r}{\sqrt{V_r}},N(0,1)\right)
\le\min\{1,C\eta_r^{1/5}\}
$$

and

$$
d_K\!\left(\frac{Y_r-\mathbb EY_r}{\sqrt{V_r}},N(0,1)\right)
\le\min\{1,C[\eta_r\log(e+\eta_r^{-1})]^{1/4}\}.
$$

The sufficient growing-threshold CLT condition improves to
$\sqrt r\,W^2/V_r\to0$. The actual-variance Bernstein MGF bound has
linear scale $2W\sqrt{\ell_r}=O(Wr^{1/4})$.

The exact constant $\ell_r$ is sharp for the deterministic two-rate
envelope inequality

$$
\int_0^\infty g_r(au)g_r(bu)\,du
\le\frac{\ell_r}{\max(a,b)},\qquad
g_r(x)=\max_{0\le k<r}\Pr(\operatorname{Poi}(x)=k).
$$

Indeed, take $a=1$ and $b\downarrow0$. This proves no necessity of a
threshold factor in a Berry–Esseen bound, no optimality of the stochastic
bracket variance, and no necessity of the displayed CLT condition.

### Why the first-discovery proof does not settle higher coverage

At threshold two, even a single species affected by independent clocks
of intensities $a,b$ has

$$
\Pr\{\operatorname{Poi}(a)+\operatorname{Poi}(b)<2\}
=e^{-(a+b)}(1+a+b),
$$

which differs from $e^{-a}(1+a)e^{-b}(1+b)$ by the extra $ab$ term in the
product. Replacing each void probability in the old hard-core transform
with a truncated Poisson factor is therefore incorrect. An expanded
representation imposes capacity constraints on total incident
multiplicities; the occupation-ratio invariant in O has not been extended
to those constraints.

A real MGF estimate alone cannot establish the missing zero-free disk:
a symmetric variable taking values $\pm\sqrt V$ is sub-Gaussian with
variance parameter $V$, yet its MGF has zeros at
$\pm i\pi/(2\sqrt V)$. This is a warning about that analytic implication,
not a counterexample in the bounded-batch model. A suitable zero-free or
variance-linear cumulant theorem for thresholds remains a concrete route
to investigate.

### A rigorous obstruction to a particular Stein-kernel route

Already for $B=1$, take $N\sim\operatorname{Poi}(r)$ and
$F=1_{\{N\ge r\}}$, and write $v_r=\operatorname{Var}F$ and
$p_r=\Pr(N=r-1)$. The add-one difference is

$$
DF(N)=F(N+1)-F(N)=1_{\{N=r-1\}}.
$$

If a proposed add-one Stein mass has the form $T=DF(N)H$ and
$\mathbb ET=v_r$, Cauchy–Schwarz forces

$$
\operatorname{Var}T\ge v_r^2(p_r^{-1}-1).
$$

For the usual conditional Mehler mass, $T$ is a function of $N$ and is
constant on its sole nonzero state, so equality holds. Since
$p_r\sim(2\pi r)^{-1/2}$ and $v_r\to1/4$,

$$
\frac{\operatorname{Var}T}{v_r}
\ge v_r(p_r^{-1}-1)
\sim\frac{\sqrt{2\pi r}}4.
$$

Thus a threshold-uniform estimate $\operatorname{Var}T\le CB^2v_r$ is
false for this add-one mass. This is **a method obstruction**: $F$ itself
is an ordinary Bernoulli variable, and sums of independent singleton
threshold indicators have the standard threshold-uniform
$C/\sqrt{\operatorname{Var}(\text{sum})}$ Berry–Esseen bound. A coupling
that replaces the indicator directly can avoid the artificial
critical-count support restriction. No counterexample to a general
threshold-uniform $CW/\sqrt{V_r}$ bound is supplied.

## 5. Signed, local, joint, and functional questions

**Signed weights and exact frequencies.** F proves variance upper bounds
and read-$B$ exponential comparisons for suitable signed or bounded
occupancy functionals. These do not yield the positive actual-variance
bracket representation used in O. Under a bound on total absolute batch
weight, the zero-free activity estimates can still be bounded using
$|h_i|$; the missing step is variance-linear control in terms of the
possibly much smaller variance of the signed sum. In particular,
$K_r-K_{r+1}$ counts species of frequency exactly $r$, but separate CLTs
for $K_r$ and $K_{r+1}$ do not prove a CLT for their difference. A signed
fluctuation theorem or a joint theorem with cancellation control is needed.

**Local limits.** A distribution supported on even integers cannot satisfy
a span-one local Gaussian approximation at odd integers. Disjoint pairs
provide this obstruction inside the model. A general result must account
for the relevant lattice or another smoothing convention and control
characteristic functions away from zero, including approximate
periodicity. The small zero-free disk in O and the endpoint-atom estimates
in Q do not supply that analysis. [Hwang–Janson][HJ] give local-limit
results for the singleton occupancy benchmark, not for arbitrary
overlapping batches.

**Joint occupancy spectra.** A vector CLT requires control of every
Cramér–Wold linear combination, including signed combinations with small
residual variance. Positive scalar results alone do not provide that
control. Mixed cumulants or a suitable multivariate coupling must be
established with the correct covariance normalization.

**Functional limits.** F supplies temporal covariance formulas and bounds
such as $\operatorname{Var}(X(t)-X(s))\le B[m(t)-m(s)]$. A process limit
still requires finite-dimensional Gaussianity under the chosen rescaling
and a tightness argument in a stated function space. Covariance
convergence and an increment second-moment bound alone are insufficient.

**Statistical calibration and random laws.** Conditional discovery
results use oracle means and variances. Unknown-law prediction intervals
need estimators, bias control, variance consistency, and appropriate
joint limits. Conditioning on a random sampling law does not by itself
give an unconditional Gaussian limit with deterministic normalization;
mixtures may retain random means or variances. These are separate
questions from the proved deterministic-law results.

## 6. Next proof obligations

1. **Fixed-sample approximation at the residual scale.** Use the now
   uniform correction to target $\tau_{h,n}^2$, and prove a residual or
   conditional approximation strong enough to recover the law at
   $N(n)=n$.
2. **Higher-coverage analytic or coupling control.** Establish a zero-free
   or cumulant bound for the capacity-constrained representation, or
   build a coupling whose mass avoids the singleton critical-count
   obstruction. Any proposed threshold-uniform estimate should pass
   the explicit $B=1$ calculation above.
3. **Signed variance and mixed cumulants.** Control cancellation before
   claiming exact-frequency, joint-spectrum, or time-indexed Gaussian
   limits.
4. **Arithmetic and process estimates.** Add the necessary off-origin
   Fourier control for local limits and the correct tightness estimates
   for a stated functional scaling.

Finite numerical checks can detect algebraic or implementation errors
in these proposals. They cannot prove uniformity over all laws,
countable-space limits, zero-free regions, or asymptotic theorems. The
repository reports record the checked finite cases and their tolerances;
the manuscripts remain the source of the mathematical claims.

## External references used to delimit the claims

- **Döring, Jansen, and Schubert.** *The method of cumulants for the normal
  approximation*. Probability Surveys **19** (2022), 185–270.
  [DOI](https://doi.org/10.1214/22-PS7) · [arXiv][DJS].
  Theorem 2.1 and Section 6 supply the classical Cramér–Petrov transfer.
- **Hwang and Janson.** *Local limit theorems for finite and infinite urn
  models*. Annals of Probability **36** (2008), 992–1022.
  [DOI](https://doi.org/10.1214/07-AOP350) · [arXiv][HJ].
  Theorems 2.1 and 2.3 and Corollary 2.5 provide singleton benchmarks.
  A local error bound must not be summed over an uncontrolled support
  and silently reinterpreted as the desired general Kolmogorov theorem.

[DJS]: https://arxiv.org/abs/2102.01459
[HJ]: https://arxiv.org/abs/math/0604397
