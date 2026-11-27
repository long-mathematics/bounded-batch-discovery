# Proof audit of the quantitative extensions — November 26, 2026

## Assessment

The mathematical review supports the three advances assembled in [Further quantitative bounds for bounded-batch discovery](../papers/quantitative_extensions.tex): relative Gaussian tails with Cramér corrections; uniform fixed-sample mean and variance corrections; and a sharper multiple-coverage bracket envelope. No unresolved mathematical gap was identified in the three reviewed arguments. Final integration, build, numerical-report, and rendering checks are recorded separately below.

This is an internal research proof audit. It is not independent peer review, formal verification, or an exhaustive novelty or priority assessment. Classical cumulant-transfer and smoothing results are attributed to their sources. The numerical examples supplement the proof review; they do not prove the infinite-species statements or universal inequalities.

The review preserves the distinction between established results and the remaining distributional problems. In particular, the fixed-sample correction theorem does **not** prove a universal fixed-sample Berry–Esseen bound, and the threshold improvement does **not** prove the first-discovery optimal rate for higher coverage.

## Model assumptions and proof dependencies

The species space is countable. Batches are independent and identically distributed random finite subsets with a deterministic law; an independent rate-one Poisson process supplies arrivals for the Poissonized statements. Empty batches are allowed. The weighted results assume deterministic nonnegative species weights whose total in every supported batch is at most a finite positive bound $W$. Individual batches remain finite, but their cardinalities need not have a uniform bound. The unweighted specialization has $W=B$.

The domination $X_h(t)\le WN(t)$, and its multiple-coverage analogue, justify all moment and exponential-moment operations. Distributional normalizations require positive variance. A zero variance is handled separately as a constant-variable case.

The relative-tail argument imports the analytic and cumulant theorem from [Optimal Berry–Esseen bounds for bounded-batch discovery](../papers/optimal_berry_esseen.tex). The fixed-sample correction and threshold arguments use the independent batch-path representation and product-space bounds from [the foundation manuscript](../papers/bounded_batch_discovery.tex). The relevant proof dependencies were rechecked, rather than inferred from the titles or numerical tests.

All three results permit changing laws, times, and weights in their stated uniform inequalities. Fixed-law asymptotic assertions elsewhere in the repository are not silently promoted to uniform changing-law statements.

## 1. Relative tails and Cramér corrections

### Analytic input and classical provenance

For $S=(X_h-\mathbb EX_h)/\sqrt V$ and $\delta=W/\sqrt V$, the imported inputs are

$$
|\kappa_k(S)|\le2k!(4\delta)^{k-2}\quad(k\ge3)
$$

and a zero-free MGF disk of radius $1/(2\delta)$. The countable-species and arbitrary-overlap extensions are already included in that input theorem.

The connection with classical Cramér–Petrov theory is correctly stated. Condition (S) and Theorem 2.1 of [Döring–Jansen–Schubert, arXiv version 2](https://arxiv.org/pdf/2102.01459v2), were checked. The conversion $\Delta=(48\delta)^{-1}$ follows from $2k(k-1)\le12^{k-2}$. Their theorem supplies the classical full correction with relative error of order $(1+x)/\Delta$. The present argument provides a direct exponential-tilt proof of its displayed constants and endpoint conventions; it does not claim the general cumulant-transfer principle as new.

### Tilt derivatives and constants

Termwise differentiation of the cumulant series gives the asserted third- and fourth-derivative bounds. On $|h|\le(1024\delta)^{-1}$, integration gives $19/20\le K''(h)\le21/20$. This proves existence and uniqueness of the saddlepoint $K'(h)=x$ for $0\le x\le(2048\delta)^{-1}$, along with the comparisons among $h$, $a=h\sqrt{K''(h)}$, and $x$.

The explicit arithmetic was checked:

| Quantity | Value | Bound used |
|---|---:|---:|
| $48(256/255)^4$ | 48.757382 | 50 |
| $768(256/255)^5$ | 783.177396 | 800 |
| Complex tilted-segment radius coefficient | 0.020706915 | $1/48=0.020833333$ |
| Complex third derivative after variance scaling | 73.419272 | 78 |
| Tilted Kolmogorov constant | 187.817974 | 200 |
| $800(20/19)^4+7500(20/19)^5$ | 10674.855892 | 11000 |

The Fourier cutoff $T=(52\delta)^{-1}$ makes the tilted Taylor remainder at most $s^2/4$. The smoothing inequality has integral coefficient $1/\pi$, integration over $[-T,T]$, and normal error $24/(\pi\sqrt{2\pi}T)$, matching [Lebowitz–Pittel–Ruelle–Speer, equation (2.16), arXiv version 2](https://arxiv.org/pdf/1408.4153v2).

### Normal integral and strict/non-strict tails

The lower bound for $J(a)=e^{a^2/2}\overline\Phi(a)$ is valid because

$$
\frac{1-e^{-1}}{\sqrt{2\pi}}=0.252179617\ldots>\frac14.
$$

For the tilted half-normal mean $m(a)$, the identities $- (\log J)'(a)=m(a)$, $m'(a)=-\operatorname{Var}_a(Y)$, and $1=a m(a)+\mathbb E_aY^2$ prove

$$
0\le-(\log J)'(a)\le\frac2{1+a},
$$

including the endpoint $a=0$.

The Laplace-transform identity uses the chosen strict or non-strict probability at zero outside the integral and $\mathbb P(U>y)$ inside it. It therefore remains correct in the presence of atoms. Kolmogorov error $\eta$ gives Laplace error at most $2\eta$ and bounds every atom by $2\eta$. Undoing the tilt proves the saddlepoint bound with constant 1800 for both signs and both endpoint conventions, without a density or lattice-span hypothesis.

The comparison of $J(a)$ and $J(x)$ is valid throughout the main range: $e^{180/2048}<1.1$. The full-correction constant 2500 follows. Under $\delta\le1/16384$ and $x\le(16384\delta)^{-1}$, its relative error is less than $1/2$, justifying the logarithmic estimates and the separate relative endpoint-atom bound.

### Legendre correction and sharp windows

The Legendre derivatives were checked directly:

$$
\Lambda'(x)=x-h(x),\qquad
\Lambda'''(x)=\frac{K'''(h)}{v^3},\qquad
\Lambda''''(x)=\frac{K''''(h)}{v^4}-\frac{3K'''(h)^2}{v^5}.
$$

They give $|\Lambda(x)|\le10\delta x^3$ and a cubic Taylor remainder at most $500\delta^2x^4$. The constants 5000 and 12000 in the logarithmic bounds then follow as stated. Replacing $S$ by $-S$ changes the cubic coefficient's sign and preserves the absolute estimates.

The resulting uniform windows are $x=o(\delta^{-1/3})$ for Gaussian relative tails, $x=o(\delta^{-1/2})$ after the cubic correction, and $x=o(\delta^{-1})$ with the full correction. The first two exponents are sharp for this model class. In the Poissonized singleton example, $X_m\sim\operatorname{Bin}(m,1/3)$, and expansion of its exact Legendre transform gives cubic coefficient $1/18$ and quartic coefficient $-1/36$. The boundary limits are respectively $e^{s^3/18}$ and $e^{-s^4/36}$, for strict and non-strict tails alike. No numerical constant or full-correction window is claimed optimal.

## 2. Uniform fixed-sample mean and variance corrections

For ordinary singleton occupancy, Hwang and Janson already prove uniform $O(1)$ mean and residual-variance corrections in Theorem 2.3 of [Local limit theorems for finite and infinite urn models](https://arxiv.org/pdf/math/0604397), *The Annals of Probability* **36** (2008), 992–1022, DOI 10.1214/07-AOP350. Their Theorem 2.1 gives a uniform local approximation and Corollary 2.5 gives a changing-law CLT for that singleton setting. These primary statements were checked and are explicitly credited in Section 7; the new theorem extends the moment corrections to overlapping batches and nonnegative weights. No extension of their local theorem is claimed.

The weighted theorem holds for every integer $n\ge1$, including deterministic or empty batches and endpoint incidence probabilities. With $f_h=m_h'$, it proves

$$
0\le b_{h,n}-m_h(n)\le2W,
\qquad
\left|s_{h,n}^2-\tau_{h,n}^2\right|
\le(4+2e^{-2})W^2,
$$

where

$$
\tau_{h,n}^2=V_h(n)-nf_h(n)^2.
$$

The proof was also independently read after its derivation. The global bound $|\delta_n'|\le2$, for $\delta_n(x)=e^{-nx}-(1-x)^n$, is valid on the full interval $[0,1]$. The covariance split includes diagonal terms and never uses $(1-a-b)^n$ when $a+b>1$.

For $\mu_1=\mathbb E\sum_{i\in A_1}h_i$ and $\mu_2=\mathbb E(\sum_{i\in A_1}h_i)^2$, the signed discrepancy is bounded by $2(\mu_2+\mu_1^2)$. The remaining exponential-series error is nonnegative and at most $2e^{-2}\mu_1^2$. This proves the sharper moment-dependent constants and asymmetric bounds in the theorem.

The residual variance is nonnegative and has both interpretations claimed:

$$
\tau_{h,n}^2
=\operatorname{Var}(X_h(n)-f_h(n)N(n))
=\min_a\operatorname{Var}(X_h(n)-aN(n)),
$$

and

$$
\tau_{h,n}^2
=n\operatorname{Var}\!\left(\sum_{i\in A_1}h_i e^{-n\lambda_i}\right)
+\sum_{i,j}h_i h_j e^{-n(\lambda_i+\lambda_j)}
 (e^{n\lambda_{ij}}-1-n\lambda_{ij}).
$$

The projection identity uses $\operatorname{Cov}(X_h(n),N(n))=nf_h(n)$. The second identity is a nonnegative decomposition, not an assumption about signs of fixed-sample species covariances.

Finite-species exhaustion gives $L^2$ convergence under the bounds $D_{h,n}\le nW$ and $X_h(n)\le WN(n)$; all discrepancy series are absolutely summable. The moment-only version requires a finite second moment of total batch weight and uses compound-Poisson domination instead. It concerns the correction identities, not a new unbounded-batch CLT.

### Sharp orders and the remaining distributional gap

For uniform sampling of $n$ singleton types at $n$ draws,

$$
b_n-m(n)\longrightarrow\frac1{2e},\qquad
s_n^2-\tau_n^2\longrightarrow\frac{3e^{-2}-e^{-1}}2>0.
$$

The expansions and diagonal terms were checked. Disjoint blocks or singleton weights scale these errors by $B,B^2$ or $W,W^2$. Thus the uniform correction orders cannot be reduced to $o(W)$ and $o(W^2)$; the numerical upper-bound constants are not asserted sharp.

A separate example shows why the Poissonized variance cannot simply replace the fixed-sample variance. For $\lceil n^{3/2}\rceil$ equally probable singleton types,

$$
V(n)\sim n,\qquad \tau_n^2\sim s_n^2\sim\frac{\sqrt n}{2}.
$$

Both variances diverge, but their ratio is unbounded. The new correction theorem supplies the correct residual scale uniformly over changing laws. It does not control the conditional distribution given $N(n)=n$, prove a normal approximation for $X_h(n)-f_h(n)N(n)$, or establish a uniform fixed-sample $CW/s_{h,n}$ Berry–Esseen bound. Those distributional steps remain unresolved here. The example is an obstruction to a transfer shortcut, not a counterexample to that conjectured rate.

## 3. Multiple-coverage envelope and inherited bounds

For integer heterogeneous thresholds $r_i\in\{1,\ldots,r\}$, the species count is $N_i(t)=\sum_{C\ni i}P_C(t)$. Both definitions are explicit in the integrated statement, as required by the integer Poisson-mass index in the jump formula. The weighted batch hypothesis is unchanged.

For the mode envelope $g_r(x)=\max_{0\le k<r}e^{-x}x^k/k!$, the exact integral is

$$
\ell_r=\int_0^\infty g_r(x)\,dx
=1+\sum_{k=1}^{r-1}\frac{e^{-k}k^k}{k!}.
$$

The proof follows by integrating $g_r-g_{r-1}=-p'_{r-1}$ on $[r-1,\infty)$. Stirling's bounds give $\ell_r\le1+\sqrt{2(r-1)/\pi}$ and $\ell_r\sim\sqrt{2r/\pi}$, with $\ell_1=1$.

Keeping the larger-rate factor and bounding the other by one gives

$$
\int_0^\infty g_r(au)g_r(bu)\,du
\le\frac{\ell_r}{\max(a,b)}.
$$

The constant is sharp for this deterministic envelope inequality: take $a=1$, $b\downarrow0$, and apply dominated convergence. This is not a lower bound on a stochastic approximation rate or proof that threshold dependence is necessary.

The terminal martingale's batch jumps lie in $[0,W]$. Its bracket has a positive pair decomposition with coefficients

$$
a_{ij}=\ell_r w_iw_j\lambda_{ij}/\max(\lambda_i,\lambda_j).
$$

The pair functions lie in $[0,1]$ and depend on independent batch paths meeting the pair. The weighted incidence identity bounds each coordinate load by $2\ell_rW^2$. Zero rates and zero coefficients are handled before division. The bracket mean is the actual variance; the bound $Q\le W^2t$ and nonnegative partial sums justify all countable weighted passages.

Putting $A=W\sqrt{\ell_r}$ gives exactly the generic foundation inputs: jumps at most $A$, bracket variance at most $2A^2V_r$, and the bracket exponential bound with load $2A^2$. The fifth-root and logarithmic fourth-root Kolmogorov bounds therefore hold with universal constants and parameter $\eta_r=\ell_rW^2/V_r$. Between arrivals, $|p_k'|\le1$ gives the derivative bounds needed for the deterministic-grid passage.

The same exponential-martingale/Hölder argument proves the actual-variance Bernstein bound with linear scale $2W\sqrt{\ell_r}=O(Wr^{1/4})$. The denominator identity and Chernoff parameter were checked. The sufficient growing-threshold condition improves to $\ell_rW^2/V_r\to0$, equivalently $\sqrt r\,W^2/V_r\to0$ when $r\to\infty$. Neither necessity of this condition nor an optimal higher-threshold $CW/\sqrt{V_r}$ bound is asserted.

## Numerical verification and its limits

The verification script [verify_quantitative_extensions.py](../scripts/verify_quantitative_extensions.py) and its report are maintained separately from this proof audit. The executed report is [quantitative_extensions_verification.json](../reports/quantitative_extensions_verification.json); its final case counts, script fingerprint, and maximum discrepancies are recorded below. The verification covers finite Poisson tilt and tail calculations, exact rational fixed-sample union laws and projection moments, threshold envelopes and weighted coordinate loads, explicit proof constants, and binomial boundary illustrations. The report distinguishes checks inside the theorem ranges from additional numerical illustrations outside those ranges.

During development, a separate union-state enumeration checked the uniform correction on 30 finite laws at sample sizes 1, 2, 3, 10, and 40. These 150 comparisons covered mean errors, asymmetric variance errors, nonnegativity, and the positive residual-variance identity; the largest identity discrepancy was $2.22\times10^{-15}$. These are preliminary development checks, not the totals of the final verification script. Their moments were obtained from union-state distributions rather than from the covariance formula being checked.

Finite enumeration, quadrature, symbolic algebra, and high-precision arithmetic can reveal implementation or algebraic errors. They cannot establish a universal statement over countably many species, prove the stated asymptotic ranges, or replace the analytic arguments above. A passing script or CI run is therefore recorded as computational evidence, not mathematical certification.

## Final integration and validation record

The detailed proof review was completed on the three component drafts, followed by a mathematical review of the integrated source. Section 7 preserves the audited variance argument and uses $T(A)$ for total batch weight; Section 8 explicitly defines the integer thresholds and species counts. The completed source credits the Hwang–Janson singleton benchmark. No mathematical gap was identified in these reviewed arguments.

The full `make all` run passed all four deterministic verification suites and all three manuscript validators. The earlier two manuscripts were rebuilt successfully and their extracted PDF text was checked against the previously published files. Their source, published PDF bytes, and historical validation reports are preserved unchanged in this update. The new manuscript contains **17 pages**; all pages were rendered with Poppler and inspected in contact sheets, with the main theorem, Legendre calculation, variance theorem, and threshold proof also examined at full-page resolution. No clipped equations, overlapping text, missing glyphs, broken headers, or other layout defect was observed.

The [document validation report](../reports/quantitative_extensions_validation.json) records all label and citation targets, the completed LaTeX log, PDF metadata, page headers, and text geometry. It reports zero errors: no unresolved or duplicate references, required reruns, TeX errors, missing characters, overfull boxes, or text-margin violations.

| Reviewed artifact | SHA-256 |
|---|---|
| Integrated LaTeX source | `39a0cb158bb57c9b250bffb613da1e09e53306b87cc479fa967a0ce1c1c6f03c` |
| Compiled PDF | `3afa735d658b297e4a80ff6c0761501506b27fb3128ef264194b4289b06044dc` |
| Verification script | `2c3b8c8c2782380cc7102e856ebf583c90a6adfd9ab03129a0263712d95fb1d4` |

The verification report also fingerprints the imported optimal-rate helper. Its calculations use 80 decimal digits, with a relative-or-absolute identity tolerance of $10^{-60}$. The largest reported high-precision identity discrepancy is below $6.93\times10^{-79}$. Nine comparisons with SciPy's double-precision binomial survival calculation have maximum log-tail discrepancy below $4.07\times10^{-15}$, against tolerance $10^{-10}$. These are numerical discrepancies, not interval-certified error estimates.

| Executed category | Cases |
|---|---:|
| Explicit finite models | 14 |
| Poisson model/time combinations | 27 |
| Exponential-tilt cases, including both signs | 216 |
| Exact rational fixed-sample laws | 70 |
| Integrated envelope thresholds | 6 |
| Two-rate envelope integrals | 30 |
| Weighted batch-coordinate loads | 216 |
| Scaled-beta comparisons with exact binomial PMFs and SciPy | 9 |
| Large-binomial boundary illustrations | 6 |

Two of the large-binomial cases lie inside the conservative full-correction theorem domain; the other four are explicitly marked as illustrations outside that finite range. At $m=10^{12}$ on the Gaussian boundary and $m=10^{16}$ on the cubic boundary, the observed full-correction relative errors are approximately $6.22\times10^{-5}$ and $7.13\times10^{-6}$, respectively. The corresponding stated theorem bounds are approximately $0.418$ and $0.364$. The scaled-beta calculation separates an analytic omitted-integral bound from the uncertified quadrature error estimate.

The new verification suite completed in about 6.8 seconds in the development environment. Its `-O` guard was checked: optimization mode exits without writing a report. The Makefile now includes all three PDFs and all four verification suites; the existing Actions workflow invokes those targets and uploads the resulting PDFs and reports. CI and merge records are available in the repository's pull-request and Actions history.

Mathematical proof review, finite numerical calculations, document validation, and rendered-page inspection are distinct checks. Their successful completion does not constitute formal certification, independent peer review, or a determination of novelty.
