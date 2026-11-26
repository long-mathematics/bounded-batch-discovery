# Audit of Sections 6–11

**Source reviewed:** `papers/bounded_batch_discovery.tex` at repository commit
`f458b6f43b5f16ddc3de9a0a0d06d0ae468785f9`, dated November 12, 2026.
**Review date:** November 26, 2026.

This is a mathematical reading and algebraic audit, not formal verification or
independent peer review. I read the complete proofs in Sections 6–11 and the
parts of Appendix A supporting the threshold and weighted extensions. I did not
run the archived numerical scripts or independently verify the bibliography;
those are separate parts of the repository audit.

## Required correction

### Specify an integer de-Poissonization window

At original line 979, Proposition `prop:quantdepois` says
`Let $V=V(n)>0$, $1\le h\le n$, and $\epsilon>0$.`
The expressions `b_{n+h}`, `b_{n-h}`, and the corresponding sample counts are
defined for integer sample sizes. Change this to:

```latex
Let $V=V(n)>0$, $h\in\{1,\ldots,n\}$, and $\epsilon>0$, with $b_0=0$.
```

This is a domain clarification, with no change to the inequality or its proof.
Defining `D_0=b_0=0` once in the model section is an equally good way to handle
the endpoint `h=n`.

## Proof checks and conclusions

No substantive mathematical error was found in the audited sections, conditional
on the product-space and bracket bounds proved earlier and on the cited standard
martingale and large-deviation theorems. The following are the checks supporting
that assessment; it is not a claim that all possible gaps have been excluded.

### Section 6: absolute variance correction

- The discrepancy lemma's bound is valid at the difficult endpoint ranges.
  For `x <= 1/2`, `-log(1-x)-x <= x^2`, and the derivative estimate is bounded
  by `2(4e^{-2}+e^{-1}) < 2`. For `x > 1/2`, the separate estimates
  `n 2^{-(n-1)} <= 1` and `n e^{-n/2} <= 2/e` hold for every integer `n >= 2`.
- With `u=a+b-c` and `w=a+b-ab`, both arguments remain in `[0,1]` and are at
  least `max(a,b)`. The covariance formula remains valid on the diagonal:
  `c=a=b` makes `u=a`, and `(1-w)^n=(1-a)^{2n}`.
- The sign in `s_n^2=V(n)+E_n-R_n` is correct. The summable majorant
  `2(c+ab)` is available exactly under the stated finite second moment.
- The remainder after removing `n f(n)^2` is bounded by
  `1/2 [sum_i lambda_i (n lambda_i) exp(-n lambda_i/2)]^2`.
  Dominated convergence applies with majorant `(2/e) lambda_i`; hence this is
  an absolute `o(1)`, rather than merely an error small relative to the variance.
- The first-moment mean comparison actually needs only `E|A_1| < infinity`,
  which is sufficient for its subsequent use in the square-root-window lemma.

### Section 7: variance equivalence and fixed-sample transfer

- A finite set of positive rates contributes negligibly to `f(t)` because a
  smaller positive rate exists outside any prescribed finite set. Cauchy–Schwarz
  on the remaining tail proves `f(t)^2/f(2t) -> 0`.
- `v_0(t)=integral_t^{2t} f(u) du >= t f(2t)` gives the claimed negligible
  correction `t f(t)^2=o(v_0(t))`. Nonnegative Poissonized covariances give
  `v_0 <= V`, also in the finite-second-moment setting.
- The split at `n^{-2/3}` works for a window of width `ceil(C sqrt(n))`:
  the small-rate multiplicative error tends to one and the large-rate part is
  smaller than every inverse power of `n`. The use of the mean comparison is
  justified by its first-moment proof, despite the broader theorem containing
  it being stated under a second moment.
- The coupling proof correctly takes `n -> infinity` for each fixed `C`, and
  only then `C -> infinity`. It does not require a uniform rate in the
  discovery-window estimate.
- The normal-density term in the quantitative transfer inequality has the
  correct coefficient and the coupling-error terms have the correct scaling.
  The integer-window clarification above is needed.
- Removing one batch yields a self-bounding function after division by `B`:
  every deletion influence lies in `[0,1]`, and their sum is at most `D_n/B`.
  The resulting bounds use the mean, with the manuscript correctly avoiding
  a claim that a fixed-sample actual-variance bound follows.

### Section 8: unbounded batches

- The truncation argument controls both the variance and the square of the
  mean of the deleted compound-Poisson incidence mass. In particular,
  `E R_l^2 / V_l <= u_l/V_l + (u_l/V_l)^2/epsilon^2` tends to zero for each
  fixed `epsilon`, which is enough for variance equivalence and centered
  `L^2` coupling.
- The retained integer batch bound is asymptotically at most
  `epsilon sqrt(V_l)`; the fourth-moment Berry–Esseen estimate therefore
  gives `C epsilon^{2/5}`. Sending the smoothing tolerance and then
  `epsilon` to zero is valid.
- In the counterexample, the later-block variance is asymptotic to
  `t_K c 2^{-K}=s_K^2`, whereas its mean divided by `s_K` is
  `O(2^{-2K})` and the probability that any later block appears tends to zero.
  Consequently the claimed two-point weak limit has variance strictly below
  one; the explanation in terms of nonuniform integrability is correct.
- The lower bound valid for every `t_K <= t < t_{K+1}` establishes divergence
  along all times, rather than only along the chosen subsequence. Replacing
  `t_K` by its integer part preserves the displayed asymptotics, so the
  fixed-sample transfer is justified.

### Section 9: actual-variance concentration and moderate deviations

- The exponential-martingale domination and Hölder exponents are correct for
  either sign of `theta`. Nonnegativity of the compensator gives the useful
  bound `Z(theta) <= exp(theta M)` also when `theta < 0`.
- With `a=B|theta|`, `q=1+a`, the denominator identity is exactly
  `1-qa/3-q^2 a/2-q(1-2a)=a/6+2a^2/3-a^3/2 >= 0` for `0<a<1/2`.
  It also verifies that the argument of the second geometric-series bound is
  below two. The final MGF inequality therefore has the actual variance as its
  quadratic coefficient.
- The Chernoff choice `|theta|=x/(V+2Bx)` gives the stated denominator
  `2(V+2Bx)`. If `V=0`, the count is constant and the tail assertion is trivial;
  an optional sentence can dispose of that case before choosing `theta`.
- The moderate-deviation proof does control the bad bracket event under the
  exponential tilt. The Cauchy–Schwarz bound is superexponentially small because
  `V_k/B_k^2 >> a_k^2`. On the good event the compensator is uniformly
  `theta_k^2 Q_k/2 (1+o(1))`, which yields the limiting log MGF.
- The theorem states an LDP and logarithmic tail limits, and correctly does
  not promote them to relative Gaussian-tail approximations or a fixed-sample
  MDP.

### Section 10: threshold counts

- Both Poisson-point envelopes are valid: the first integrates to the factor
  `2^{2r-1}/(lambda_i+lambda_j)`, and the second to
  `e^2 r/(lambda_i+lambda_j)`. For `r=1` the exact discovery envelope gives
  `d_1=1`.
- Using entire incident Poisson paths makes the support/load argument valid
  for thresholds. Using only first-arrival clocks here would be incorrect,
  and the manuscript explicitly avoids that error.
- The factor `r` enters the fourth-moment estimate through the bracket-load
  bound; the jump bound remains `B`, giving an absolute constant in the
  growing-threshold CLT condition `rB^2/V_r -> 0`.
- The logarithmic refinement uses the same completion error, which is at most
  `C(rB^2/V_r)^{2/7}`, and bracket tails on the correct scale.
- For fixed `r`, a bracket load at most `2e^2 r B^2` changes constants in the
  small-MGF estimate but leaves the stated moderate-deviation range intact.
  The manuscript does not claim a growing-`r` MDP under the fixed-`r` condition.
- Appendix A's derivative bound is valid uniformly in `r`: the absolute
  difference of two Poisson point probabilities is at most one.

### Section 11: weighted, conditional, and temporal statements

- For nonnegative weights, `X_w(t) <= W N(t)` supplies every moment even if
  the unweighted batch cardinality has no uniform bound or finite mean.
  Tonelli and the within-batch weight bound give
  `sum_j w_j lambda_ij <= W lambda_i`; the bracket-coordinate load is at
  most `2W^2`. Thus the same proofs apply for every real `W>0`.
- The phrase “with `B` replaced by `W`” refers to the listed normal-approximation,
  concentration, and MDP results, all of which use the bound `2B^2`; it must not
  be read as substituting `W` into the separate combinatorial formula
  `B(4B-3)/2`. The present theorem is clear enough in context, though this
  distinction can be made explicit in documentation.
- Conditioning on the past leaves independent future Poisson increments.
  Intersecting each batch with the now-fixed unseen set produces exactly the
  permitted bounded-batch model; the universal constant makes the conditional
  normal approximation valid for almost every history.
- The joint-absence formula across times gives the stated temporal covariance
  with `min(s,t) lambda_ij`. The increment indicators are functions of the
  incident batch clocks and satisfy the read-`B` variance bound without any
  monotonicity or association claim.

## Scope of this finding

These checks support retaining the stated results in Sections 6–11, with the
integer-window clarification. They do not independently establish the earlier
bracket inequalities, the external normal-approximation theorems, novelty or
priority, or the accuracy of archived reproducibility claims. Those items need
the corresponding core, reference, and script audits.
