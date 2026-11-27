# Fixed-sample proof audit — November 26, 2026

## Assessment and exact status

The internal review supports the proofs in [Fixed-sample zero-free bounds for bounded-batch discovery](../papers/fixed_sample_bounds.tex). Two separate proof readers checked the graph induction, limiting representation, harmonic normal-approximation argument, Hoeffding-kernel operator estimate, sector obstruction, and lattice lower bound. No mathematical gap was identified in those stated results.

The full uniform optimal fixed-sample Berry–Esseen bound **is not proved**. The new theorem has a factor `log(n+2)`. An optimal-order bound is proved when the variance is at least a fixed positive multiple of `n W²`; the constant depends on that multiple. The sufficient actual-variance MGF estimate in Proposition 7.1 remains an unproved hypothesis for the full class. The manuscript, README, and research-status note retain this distinction.

This is an internal research audit, not independent peer review, formal verification, or an exhaustive novelty assessment. The harmonic normal-approximation input is explicitly credited to Michelen–Sahasrabudhe. Deterministic finite calculations supplement the proof review and do not establish continuum nonvanishing or universal asymptotic claims.

## 1. Model and limiting operations

The main disk and logarithmic-rate theorems allow independent, nonidentically distributed random finite batches on a countable species space. Real species weights satisfy

$$
\sum_{i\in A_k}|h_i|\le W
$$

for every draw, almost surely. Thus `|D| ≤ nW`; no positive lower bound on weights, lattice mesh, number of species, or species dependency degree is assumed. Individual batches are finite, while their cardinalities may be unbounded. Empty batches and deterministic choices cause no problem. The normal-approximation statements require positive variance.

The Hoeffding results explicitly specialize to iid draws and nonnegative weights. The sector counterexample and lower-bound example are iid and unweighted. These assumptions are not silently transferred to or from the more general main theorem.

Finite species projection preserves the absolute-weight bound and each draw's independence. For a finite species set, each projected draw has a finite law automatically. The limit in the auxiliary activity parameter is a finite coefficientwise polynomial limit, locally uniform after the MGF substitution. The limiting function has value one at zero, so Hurwitz excludes an identically zero limit. The countable-species passage uses locally uniform convergence of MGFs from `|D| ≤ nW` and a second application of Hurwitz. The same domination makes the Hoeffding projection limits valid in `L²`.

## 2. Clique occupation-ratio induction

The graph has an independent species side; the other side induces a disjoint union of cliques. A right clique represents the mutually exclusive batch options for one draw. Right activities are nonnegative real numbers; left activities may be complex.

The induction is simultaneous over all induced subgraphs: nonzero partition functions, left ratios of modulus at most their activity radii, and right ratios with nonnegative real part. Every denominator used in a recursion belongs to a proper induced subgraph.

For a left vertex, the usual deletion product contains right factors `1+R_j`, each of modulus at least one. For a right vertex in clique `C`, removing the whole clique first gives

$$
R_j=\frac{w_jq_j}{1+\sum_{k\in C\setminus\{j\}}w_kq_k}.
$$

Each `q_k` telescopes through left deletions and has argument bounded by the sum of the corresponding `arcsin(r_i)`. With budget strictly less than `π/4`, both numerator and denominator lie in that cone, so the quotient has argument strictly between `−π/2` and `π/2`. The constant term one makes the denominator nonzero even when all other activities vanish. Zero right activities yield zero ratios directly.

The `π/4` budget is material. Copying the bipartite `π/2` invariant would fail because the denominator in the categorical ratio also rotates. A one-draw example in the research notes verifies that the old stronger ratio invariant need not hold. The published proof uses only the valid `π/4` form.

With right activities `M μ_k(C)`, expanding by the species subset gives exactly

$$
\mathcal Z_M(u)=\sum_{S}u_S\prod_k[1+M\Pr(A_k\cap S=\varnothing)].
$$

Division by `(1+M)^n` and passage to infinity recovers the fixed-sample absence transform. For `u_i=e^{-zh_i}−1`, the inequalities `e^a−1≤2a` and `arcsin r≤πr/2` give the angle budget when `|z|W<1/4`. The argument uses absolute weights and is therefore valid for signed weights as stated.

## 3. The logarithmic rate and its primary-source input

The primary source checked was Marcus Michelen and Julian Sahasrabudhe, *Central limit theorems and the geometry of polynomials*, **J. Eur. Math. Soc. 28** (2026), 2261–2305, [DOI 10.4171/JEMS/1530](https://doi.org/10.4171/JEMS/1530), including the [official full text](https://ems.press/content/serial-article-files/52506?nt=1).

The imported statements are Lemma 4.1 and Corollary 9.3, with the real-valued extension explained in Section 12. The definition of the normalized coefficients was checked: `a_j=κ_j/j!`, so `a_2=σ²/2`. The argument does not assume that arbitrary real weights produce a polynomial PGF.

For `Y∈[0,N]` with zero-free MGF on `|z|<1/4`, the positive harmonic function `N/4−Re K(z)` has value `N/4` at zero. Harnack gives `|Re K|≤N/2` on `|z|≤1/8`. This is the model-specific replacement for a root-count/degree estimate.

The explicit parameters are

$$
R=e^{1/16},\quad r=e^{1/32},\quad
\delta=[64\log(2N+2)]^{-1},\quad \epsilon=\delta/4.
$$

The truncated sector lies strictly inside the known harmonic region, and

$$
(r/R)^{1/\delta}\max |u|
\le\frac{N}{2(2N+2)^2}\le\frac1{32}<\frac38.
$$

This verifies Lemma 4.1 with `b=1`. The disk of radius `ε` around one lies inside the resulting smaller sector. If `σ>64/ε`, the variance term alone exceeds one in Corollary 9.3's coefficient sum; otherwise the bound follows from `d_K≤1` with a larger universal constant. Scaling `D/W` gives `log(n+2)`; signed weights are handled by translating `D/W∈[−n,n]` into `[0,2n]`.

No explicit small numerical Berry–Esseen constant is claimed. In particular, the Poissonized constant 150 is not transferred to this theorem. Michelen–Sahasrabudhe's disk sharpness result concerns general distributions, not a lower bound requiring a logarithm for the discovery class.

## 4. Canonical kernels and operator norms

The orthogonal expansion is

$$
D-\mathbb ED=\sum_{k=1}^n\sum_{|J|=k}H_k(A_J),\qquad
H_k=(-1)^{k+1}\sum_i h_iq_i^{n-k}\xi_i^{\otimes k},
$$

where `ξ_i=1{i∈A}−p_i`. Different coordinate subsets are orthogonal by conditioning on a coordinate in their symmetric difference. Thus the variance is the sum of `binom(n,k) ||H_k||²_2`.

For a canonical tensor `g`, centering in every coordinate allows `ξ_i` to be replaced by hit indicators inside its inner product. Cauchy–Schwarz supplies a factor `p_i^r`, which cancels the factor `p_i^{−r}` in `||T_r g||²`. The remaining species sum is the weight of an intersection of `r` batches, at most `W`. Therefore `||T_r||≤√W`.

The flattening factorization is

$$
H_k=\pm T_{k-r}^*\operatorname{diag}(p_i^{k/2}q_i^{n-k})T_r.
$$

The binomial probability bound then proves `sqrt(binom(n,k)) ||H_k||op≤W`. Zero-probability species contribute nothing; probabilities one give zero centered incidence. Countable limits follow by continuity of orthogonal projections and domination by `nW`.

The geometric refinement for `p_i≤1/2` uses the exact identity

$$
\binom nkp_i^kq_i^{2(n-k)}
=(1+q_i)^{-k}\Pr\{\operatorname{Bin}(n,1-q_i^2)=k\}
\le(2/3)^k.
$$

Species with `p_i>1/2` have total weight at most `2W`. For that part, the rank-one norm `||ξ_i||²=p_iq_i` yields the additional bound `2W2^(−n/2)`. The full bounds are therefore summable over all orders with any fixed polynomial weight, as stated in Corollary 5.2. The manuscript does not confuse an operator norm with the Hilbert–Schmidt norm defining the variance, or claim an unproved bound on mixed-order contractions.

## 5. Sector obstruction and sharp lower order

The sector construction takes two communities of `M≥3` species, with each draw a uniform pair in a uniformly selected community. The capacity bound

$$
\frac{\binom a2+\binom b2}{2\binom M2}
\le2^{(a+b)/M-2}
$$

has equality only for `(a,b)=(M,0),(0,M),(M,M)`. The equality characterization was checked, including zero-capacity sets and boundary values. It proves that only degrees `M` and `2M` survive the rescaled PGF limit. Direct coverage probabilities give their limiting coefficients two and one. Rouché applies at a simple nonzero root of the limiting polynomial.

For fixed `M`, individual missing species dominate and the variance is asymptotic to `2M(1−1/M)^n`, tending to zero. The roots move towards the positive PGF axis near zero, and their logarithms move far left in the MGF plane. This is consistent with the local disk and supplies no counterexample to the optimal normal-approximation rate. A tilted distribution is not assumed to be another iid batch model.

The separate lower-order example uses uniform sampling from `m` disjoint blocks at `n=m`. Its variance grows like `B²m(e^{−1}−2e^{−2})`. The lattice and Chebyshev argument gives maximum atom at least `3/[4(4s/B+2)]`. Half that jump is a lower bound against any continuous CDF, yielding `B/(16s)` for `s≥B`. This establishes the necessity of the order `B/s`, not the necessity of the logarithm in the upper bound.

## 6. Sufficient optimal criterion and the proved large-variance regime

Proposition 7.1 assumes a real centered MGF bound `A s²θ²` on `|θ|≤a/W`. On the complex circle of radius `r=min(a/2,1/8)/W`, the modulus inequality supplies the same upper bound for the real part of the analytic logarithm. Its harmonic Fourier coefficients yield

$$
|\kappa_k|/k!\le2A s^2r^{2-k}.
$$

The displayed bound is for `k≥2`. Integration review corrected the initial draft's `k≥1` wording: the analytic logarithm is centered, so its first coefficient is zero, whereas the uncentered first cumulant is the mean. No first-cumulant estimate is used in the argument.

At cutoff `|t|≤sr/[16(A+1)]`, the Taylor remainder divided by `t²` is at most `2A/[15(A+1)]<1/4`. The Gaussian-damped smoothing integral therefore gives `C(A,a) W/s`. This conditional implication is proved, but its hypothesis is not asserted for all fixed-sample models.

For `s²≥η nW²`, bounded differences does provide the hypothesis. Replacement of one batch has range at most `W` for nonnegative weights and `2W` for signed weights. Conditional Hoeffding bounds along the Doob martingale give real log MGF bounds `nW²θ²/8` and `nW²θ²/2`. Thus the optimal order holds with a constant depending on `η`. The full problem permits `s²/(nW²)→0`, and this dependence cannot be suppressed.

## 7. Build and computational verification

The [verification script](../scripts/verify_fixed_sample.py) and its [generated report](../reports/fixed_sample_verification.json) record the finite models, arithmetic precision, tolerances, case counts, and source fingerprint. Probabilities use exact rational or integer enumeration; complex values use high-precision arithmetic. Finite complex grids are explicitly not nonvanishing certificates. The two-community coefficient calculation compares independent dynamic programming and inclusion–exclusion formulas.

The manuscript is compiled from the committed LaTeX source. The [document validation report](../reports/fixed_sample_validation.json) records citation and label checks, the final compilation log, PDF text, margin geometry, and source/PDF fingerprints. The existing manuscripts and their historical validation reports remain unchanged. The full local build passed all five verification suites and all four manuscript validators. The new suite records 82,846 checks, including two signed-weight models and 18 numerical SVD cases. All ten rendered pages were visually reviewed, with the harmonic proof, summable operator estimate, and optimal-rate criterion also inspected at full page size. Repository CI is required to pass before merging the associated research pull request.
