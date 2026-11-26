# Core proof audit: Sections 2–5 and Appendix A

Audit date: 2026-11-26. Source examined: `papers/bounded_batch_discovery.tex`,
repository head `f458b6f43b5f16ddc3de9a0a0d06d0ae468785f9`, before the repository
revision. Line numbers below refer to that source. This is a proof review, not a
formal verification or independent peer review. Numerical checks were not used
as substitutes for the arguments below.

## Result

I found no substantive mathematical error in the reviewed arguments for the
product inequalities, variance refinements, bracket decomposition, normal
approximation, smoothing, Gaussian completion, or deterministic-grid passage.
One harmless degenerate-case clarification is warranted: the exponential
weighted-overlap proposition divides by its load parameter and should explicitly
require `L>0` (or declare the continuous interpretation at `L=0`). No application
in the manuscript is affected, because its load parameters are strictly positive.

The parameter range `p>=1` in the quotation of Haeusler's inequality is correct. Mourrat's Theorem 1.1 explicitly gives `p in [1,+infinity)`.
The proof here only uses fixed `p>=2` in any event.

## Precise recommended edit

At lines 253–256, replace

> Under (the overlap condition), suppose ...

with

> Under (the overlap condition), with `L>0`, suppose ...

Alternatively append that the case `L=0` is interpreted by continuity and has a
constant sum. Requiring `L>0` is simpler. Do not alter the variance-drop lemma:
its inequality is meaningful at `L=0`.

## Verification details

### Section 2: product-space inequalities (lines 175–300)

- The martingale-difference proof of weighted variance drop is valid. If a
  coordinate is outside a summand's support, independence makes the corresponding
  difference zero even after the earlier coordinates have been exposed. Weighted
  Cauchy–Schwarz therefore has precisely the displayed coefficient load, and
  orthogonality gives the claimed inequality.
- The countable nonnegative extension is justified by the uniform second-moment
  bound for partial sums and monotone convergence. For centered sums with finite
  right side, applying the same inequality to tails yields the requisite L2
  convergence.
- The Finner proof handles countably many independent coordinates correctly:
  truncate each nonnegative integrable function, condition on the first K
  coordinates, use support preservation under independence, then take K to
  infinity and the truncation level to infinity in that order. Constant functions
  with empty supports may indeed be factored out.
- For the exponential overlap bound, convexity gives
  `E exp(L theta F) <= 1+(exp(L theta)-1) E F` for either sign of theta.
  Applying `log(1+u)<=u` is legitimate because the quantity inside the logarithm
  is positive. Countable passage before centering is valid for both signs as
  described.
- For a centered indicator U, `E|U|^k <= E U^2` for all integer k>=2. This gives
  the stated two-sided mgf estimate, followed by the displayed Bennett rate.
  The variance proxy in that calculation is `B v`, consistently with the tail.

### Section 3: variance and direct concentration (lines 302–453)

- The covariance expression includes the diagonal correctly; nonnegativity
  gives the lower comparison and read-B gives the upper comparison. Singleton
  laws and disjoint B-block laws attain the constants 1 and B respectively.
- Signed weights are correctly absorbed into the individual functions; this
  argument requires no sign assumption on the weights. The manuscript expressly
  separates centered L2 convergence from existence of an uncentered sum.
- For the species-dependent refinement, applying Cauchy–Schwarz with the factor
  `|C|` is valid even when only a subset of C belongs to the finite test set J.
  The resulting integral is exactly
  `integral_0^t lambda_i exp(-2 lambda_i(t-s)) exp(-lambda_i s) ds`
  `= exp(-lambda_i t)(1-exp(-lambda_i t))`.
  No finite global batch-size moment is needed for a finite J; infinite kappa_i
  makes the bound vacuous, as the statement explains.
- Fractional Bennett uses the correct reciprocal scaling in the exponential.
  The special choice `alpha_i=1/B_i` is feasible because `B_i>=|C|` for every
  i in a supported batch C. Fatou gives the countable mgf bound without a hidden
  exponential-integrability assumption.
- The direct Bernstein denominator `2Bv_0+(2B/3)x` follows from
  `h(u)>=u^2/(2(1+u/3))`; no factor is lost.

### Section 4: terminal martingale and bracket (lines 455–584)

- Conditioning on all arrivals up to s yields the displayed Doob martingale.
  Its jumps lie in [0,B]; its compensated representation and the identity
  `E Q_t=V(t)` hold because the aggregate arrival rate is one and all integrands
  are bounded. Empty batches cause only zero jumps.
- Integrating each ordered pair gives exactly `a_ij F_ij`, where
  `a_ij=lambda_ij/(lambda_i+lambda_j)`. Each F is nonnegative and at most one.
  Its independent-coordinate support is the collection of batch clocks that
  meet either species.
- The load computation uses ordered pairs correctly. With
  `a_ii=1/2` and `sum_(j!=i) a_ij<=B-1`, its value is at most
  `2|C|(B-1/2)-|C|/2 = |C|(2B-3/2) <= B(4B-3)/2`.
  The diagonal subtraction is essential and present.
- The sharp-load construction checks out: there are B diagonal hub terms,
  `2BN(B-1)` ordered hub–leaf terms of size `1/(N+2)`, and `B(B-1)` ordered
  distinct-hub terms of size `1/(2(N+1))`. Their sum tends to the stated L_B.
  The manuscript properly limits this sharpness claim to this coefficient-load
  representation, without claiming that the actual bracket variance constant
  is optimal.
- Chernoff optimization of the bracket mgf gives both displayed tails, and
  integration of the associated Bernstein bound gives the fixed-p moment bound.

### Section 5: normal approximation (lines 586–763)

- Haeusler's moment expression is accurately quoted from Mourrat's Theorem 1.1.
  With variance normalized to one it is the sum of the p-th bracket error and
  the summed 2p-th increment moments, raised to `1/(2p+1)`. Bounded discrete-grid
  increments are not required.
- Compensator integration gives the jump-moment bound
  `E sum |Delta M|^(2p) <= B^(2p-2)V`. For eta=B²/V<=1, the three terms are
  `eta^(p/2)`, `eta^p`, and `eta^(p-1)`. Since p>=2, they are bounded by a
  constant times `eta^(p/2)`. This yields exponent `p/(4p+2)` as stated.
- The smoothing lemma does not require independence of Z and sigma. Its proof
  first bounds the mass of Z in any interval of length a on `{sigma<=a}` by
  the mass of Y in the interval enlarged by a on both sides. Independence of
  the added G is sufficient. Partitioning each half-line into half-open
  intervals controls the crossing probability by the summable Gaussian tails.
  Atoms and endpoints are covered by the one-sided-limit sentence.
- The stopped bracket is `min(Q_t,nu)` because Q is continuous and increasing.
  Appending the conditionally scaled independent Gaussian increments has mean
  zero and unconditional total variance nu. Its variance contribution is
  known at the start of the appended block. The sixth-moment contribution is
  exactly `15 m^(-2) E(nu-Q_tau)^3`, and the stopped martingale contributes at
  most `B^4 nu`. Thus the completion error has the asserted power
  `(B^4/nu²)^(1/7)`.
- For nu=(1+delta)V, scale change costs O(delta). On the good bracket event,
  added noise has standard deviation at most sqrt(2delta) after normalization
  by sqrt(V). Stopping changes the terminal variable only on the upper-tail
  event. The two bracket tails yield `exp(-c delta²/eta)`. The final bound
  `eta^(2/7)+sqrt(delta)+exp(-c delta²/eta)` is therefore justified.
  Substituting `delta=A sqrt(eta log(1/eta))` yields the fourth-root logarithmic
  rate. The completion term is smaller since `2/7>1/4`.

### Appendix A: grid passage (lines 1623–1725)

- Between arrivals the discovery density has derivative bounded by `2B²`.
  The stated uniform threshold estimate follows from
  `p_k'=p_(k-1)-p_k` and `|p_k'|<=1`; weighted discovery uses the total weight W.
- Conditional martingale isometry gives the grid predictable variance. The
  decomposition into the no-arrival event and the event of an arrival yields
  the asserted Lp density error, and the conditional-expectation contraction
  then gives Lp convergence to Q_t.
- Finite activity makes the high-power variation converge pathwise to the
  summed jump powers. The Poisson total-variation bound gives an integrable
  dominating variable for every fixed p. Grid sums equal M_t exactly, so there
  is no additional terminal-law convergence step to justify.
- The stopped-density argument correctly avoids claiming an O(mesh) crossing
  probability for the stopping time. Piecewise continuous bounded sample paths
  give left-grid convergence for Lebesgue-almost every time, and product-space
  dominated convergence plus Hölder gives the integrated Lp convergence needed.
  The stopped total variation remains bounded by the same dominating variable.
- The model-by-model grid limit does not impair uniformity of the final
  Berry–Esseen constant, because the latter comes from the uniform bracket and
  jump bounds after the grid limit.

## Primary external source checked

Jean-Christophe Mourrat, *On the rate of convergence in the martingale central
limit theorem*, author's PDF, Theorem 1.1 and the augmentation argument in the
proof of Theorem 1.5:

https://perso.ens-lyon.fr/jean-christophe.mourrat/martingale_CLT.pdf

Theorem 1.1 states p in [1,+infinity); the displayed formula matches the present
manuscript after setting its unconditional variance to one. Some abstract/search
text renders the range inconsistently, so the theorem statement is the relevant
evidence. No unverified claim about a newer theorem is needed for Sections 2–5.
