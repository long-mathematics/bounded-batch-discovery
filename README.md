# Bounded-batch species discovery

**Optimal Poissonized Berry–Esseen bounds, fixed-sample bounds with a logarithmic factor, and sharp relative-tail windows for discovery through random batches.**

Christopher D. Long — research drafts, November 2026.

[Higher coverage (PDF)](papers/higher_coverage.pdf) · [LaTeX source](papers/higher_coverage.tex)

[Indicator switching for pairs (PDF)](papers/pair_indicator_coupling.pdf) · [LaTeX source](papers/pair_indicator_coupling.tex) · [Proof audit for both papers](audit/higher_coverage_and_pairs_2026-11-27.md)

[Fixed-sample paper (PDF)](papers/fixed_sample_bounds.pdf) · [LaTeX source](papers/fixed_sample_bounds.tex) · [Proof audit](audit/fixed_sample_2026-11-26.md)

[Further quantitative bounds (PDF)](papers/quantitative_extensions.pdf) · [LaTeX source](papers/quantitative_extensions.tex) · [Proof audit](audit/quantitative_extensions_2026-11-26.md) · [Research status and open questions](extensions/research_status.md)

[Optimal-rate paper (PDF)](papers/optimal_berry_esseen.pdf) · [LaTeX source](papers/optimal_berry_esseen.tex) · [Proof audit](audit/optimal_rate_2026-11-26.md)

[Foundation paper (PDF)](papers/bounded_batch_discovery.pdf) · [LaTeX source](papers/bounded_batch_discovery.tex) · [Foundation audit](audit/2026-11-26.md) · [Original bundle](archive/2026-11-12/README.md)

This project studies the number of distinct species discovered through independent random batches. Species within a batch may be dependent. A batch may be a graph edge, a bounded hyperedge, or another finite set of observations.

The November 26 optimal-rate paper proves the original target in this repository: a uniform Poissonized Kolmogorov bound of order $B/\sqrt V$, with the explicit constant 150. The follow-up paper proves relative Gaussian tails and a cubic correction with sharp ranges, uniform fixed-sample mean and variance corrections, and sharper bounds for higher coverage. The new fixed-sample paper proves a zero-free MGF disk and the uniform rate $CW\log(n+2)/s_{h,n}$, including signed weights and independent batches with different laws. The November 27 papers prove a higher-coverage rate of $15W\sqrt{\ell_r}/\sqrt V$ for arbitrary overlap and a threshold-uniform $CW/\sqrt V$ for star forests. Removing the higher-coverage threshold factor on arbitrary pair graphs and removing the fixed-sample logarithm remain separate open questions. Internal mathematical audits checked the arguments and their imported estimates; the papers remain research drafts, without independent peer review or formal verification. The linked audit records describe the checks and their limits.

## Main results

For i.i.d. batches $A_1,A_2,\ldots$ on a countable species space, with a deterministic law and $|A_k|\le B$, write

$$
D_n=\left|\bigcup_{k=1}^{n}A_k\right|,\qquad
X(t)=D_{N(t)},\qquad V(t)=\operatorname{Var}X(t),
$$

where $N$ is an independent rate-one Poisson process. Set $\lambda_i=\Pr(i\in A_1)$ and $v_0(t)=\sum_i e^{-t\lambda_i}(1-e^{-t\lambda_i})$.

| Result | Scope and manuscript location |
|---|---|
| Optimal-order Berry–Esseen bound $\min\{1,150B/\sqrt V\}$ | All Poissonized bounded-batch laws and times with $V>0$; no dependency-degree bound. Companion Theorem 1.1. |
| Weighted bound $\min\{1,150W/\sqrt{V_w}\}$ and variance-linear cumulants | Nonnegative weights with total weight at most $W$ per batch; cardinalities may be unbounded. Companion Theorem 1.2. |
| Relative Gaussian tails for $x=o(\delta^{-1/3})$; cubic correction for $x=o(\delta^{-1/2})$ | Poissonized first discovery, $\delta=W/\sqrt V$; both ranges are sharp. Full Cramér correction has relative error $O(\delta(1+x))$ for $x\le c/\delta$. Follow-up Sections 1–6. |
| Uniform $O(W)$ mean and $O(W^2)$ variance corrections | Fixed sample size, including changing laws; residual variance $V_w(n)-nf_w(n)^2$. The error orders are sharp. Follow-up Section 7. |
| Fixed-sample zero-free MGF on $\lvert z\rvert<1/(4W)$ and rate $CW\log(n+2)/s_{h,n}$ | Independent batches with possibly different laws; countable species, arbitrary overlap, signed weights with absolute batch total at most $W$. Fixed-sample paper, Theorems 1.1–1.2. |
| Uniform, summable bounds for Hoeffding-kernel flattenings | Iid nonnegative-weight batches. Fixed-sample paper, Proposition 5.1 and Corollary 5.2. Mixed-contraction control needed for the optimal rate remains open. |
| Optimal fixed-sample order when $s_{h,n}^2\ge\eta nW^2$ | Constant $C_\eta$ depends on the positive variance ratio. Fixed-sample paper, Corollary 7.2. |
| No universal zero-free PGF sector about the positive real axis | A two-community construction with $B=2$. Fixed-sample paper, Theorem 6.1. This does not refute the optimal Berry–Esseen conjecture. |
| Threshold bracket load $2\ell_rW^2$, where $\ell_r\sim\sqrt{2r/\pi}$ | Nonnegative weights and heterogeneous integer thresholds at most $r$; improves generic normal approximation and actual-variance concentration. Follow-up Section 8. |
| Higher-coverage Kolmogorov error $\le 15W\sqrt{\ell_r}/\sqrt V$ | Nonnegative weights, arbitrary overlap, integer thresholds at most $r$. The variance exponent is sharp for each fixed threshold. |
| Threshold-uniform error $\le CW/\sqrt V$ on star forests | Pair and singleton batches, heterogeneous thresholds and component sizes. On arbitrary pair graphs the aggregate squared switching cost is at most $3WV$ and $\lvert\kappa_3(Y)\rvert\le4WV$. |
| Variance-normalized central limit theorem whenever $\operatorname{Var}(D_n)\to\infty$ | Fixed law; no connectivity, regular variation, or within-batch independence assumption. Foundation Theorem 1.1 and Section 7. |
| Sharp variance comparison $v_0(t)\le V(t)\le Bv_0(t)$ | Poissonized discoveries. Foundation Section 3. |
| Actual-variance Bernstein concentration and Gaussian moderate deviations | Poissonized counts. Foundation Section 9. |
| Absolute expansion $\operatorname{Var}(D_n)=V(n)-n m'(n)^2+o(1)$ | Fixed laws with $\mathbb E|A_1|^2<\infty$, where $m(t)=\mathbb E X(t)$. Foundation Section 6. |
| Lindeberg extension and a finite-second-moment counterexample | Unbounded batches; finite second moment alone does not imply a CLT. Foundation Section 8. |
| Threshold, conditional, and temporal results | Each with its stated assumptions. Foundation Sections 10–11. The higher-coverage papers give the quantitative results above. |

The new bound is

$$
d_K\!\left(\frac{X(t)-\mathbb E X(t)}{\sqrt{V(t)}},Z\right)
\le \min\!\left\{1,\frac{150B}{\sqrt{V(t)}}\right\},
\qquad Z\sim\mathcal N(0,1).
$$

The constant is uniform in the law, time, number of species, and overlap pattern, so it applies to sequences of changing laws and batch bounds. The order is sharp: disjoint batches give $X=B\operatorname{Bin}(m,1/2)$, whose largest atom has order $B/\sqrt V$. The numerical constant is not claimed optimal; the proof actually gives $36/\sqrt\pi+864/(\pi\sqrt{2\pi})<131$.

The foundation paper supplies the actual-variance exponential-moment inequality used in the new proof, as well as

$$
\Pr\bigl(|X(t)-\mathbb E X(t)|\ge x\bigr)
\le 2\exp\!\left(-\frac{x^2}{2(V(t)+2Bx)}\right),\qquad x>0.
$$

For $B=2$, the fixed-law theorem addresses the vertex-count question in Janson's countable-edge model without the connectivity hypothesis used in Eriksson's connected-regime result. The [reference audit](audit/references_2026-11-26.md) checks the cited problem and theorem statements. Diffuse dust mechanisms and arbitrary mixtures over random laws are outside this claim.

## How the optimal rate is proved

The new paper combines two properties that are uniform in the number of species:

1. The absence polynomial is exactly a normalized independence polynomial of the bipartite species–batch incidence graph. A simultaneous occupation-ratio induction proves that the discovery MGF has no zeros for $|z|<1/(2B)$. Only batch sizes enter; species may have arbitrarily large dependency degree.
2. The foundation paper proves $\log\mathbb E e^{\theta(X-\mathbb EX)}\le V\theta^2/[2(1-2B|\theta|)]$. This controls the real part of the analytic log MGF on the smaller disk $|z|\le1/(4B)$.

Taking harmonic Fourier coefficients gives

$$
|\kappa_k(X)|\le 2k!V(4B)^{k-2},\qquad k\ge3.
$$

Esseen's smoothing inequality then gives the stated Kolmogorov rate. The same argument uses the total batch weight $W$ in the weighted model. Zero-freeness alone is insufficient for the variance-linear cumulant estimate; the actual-variance moment bound is essential.

[Berry–Esseen bounds beyond sparse dependency graphs](extensions/berry_esseen.md) retains the original comparison with Maximilian Janisch's dependency-graph manuscript and records how the new proof resolves the target. The new proof does not depend on that external theorem or its Lean formalization. The foundation paper retains its earlier valid, weaker quantitative bounds; the companion paper supplies the stronger first-discovery result.

## Further results and what remains open

The [follow-up paper](papers/quantitative_extensions.pdf) develops three consequences and extensions. The [research-status note](extensions/research_status.md) gives precise hypotheses, remaining proof obligations, and unsuccessful routes.

### A fixed-sample rate, with the logarithm still to remove

For independent draws $A_1,\ldots,A_n$ satisfying $\sum_{i\in A_k}|h_i|\le W$, set $s_{h,n}^2=\operatorname{Var}D_{h,n}>0$. The [fixed-sample paper](papers/fixed_sample_bounds.pdf) proves

$$
\mathbb Ee^{zD_{h,n}}\ne0\quad (|z|<1/(4W)),\qquad
d_K\!\left(\frac{D_{h,n}-\mathbb ED_{h,n}}{s_{h,n}},Z\right)
\le\min\left\{1,\frac{CW\log(n+2)}{s_{h,n}}\right\}.
$$

Each draw is represented by a clique of mutually exclusive batch choices. A new occupation-ratio induction proves the disk; Michelen–Sahasrabudhe's harmonic normal-approximation theorem, with a bounded-support Harnack estimate, gives the rate. No lattice mesh, species-degree bound, or identical-law assumption is required.

For iid nonnegative weights, every flattening of the canonical Hoeffding kernel obeys $\sqrt{\binom nk}\|H_k\|_{\mathrm{op}}\le W$. A refinement makes these bounds summable over the orders, even with any fixed polynomial weight. This does not yet control their mixed contractions. The paper also proves that a global zero-free PGF sector fails even for $B=2$. Its roots approach the positive axis near zero, while the variance tends to zero; this is an obstruction to that analytic shortcut, not to the desired rate.

The uniform target $CW/s_{h,n}$ remains open here. The order is necessary by a disjoint-block lattice example. A local bound $\log\mathbb E e^{\theta(D-\mathbb ED)}\le A s_{h,n}^2\theta^2$ for $|\theta|\le a/W$, with universal $A,a>0$, would remove the logarithm using the proved disk. Bounded differences supplies it when $s_{h,n}^2\ge\eta nW^2$, proving the optimal order with a constant $C_\eta$. The full variance-sensitive estimate, without that restriction, has not been established.

### Relative tails, including lattice endpoints

For $S=(X_w-\mathbb EX_w)/\sqrt V$ and $\delta=W/\sqrt V$, let $K(h)=\log\mathbb Ee^{hS}$. Given $x\ge0$, choose $h$ with $K'(h)=x$ and put $\Lambda(x)=x^2/2-hx+K(h)$. The proved cumulant bounds yield

$$
\Pr(S\ge x)=\overline\Phi(x)e^{\Lambda(x)}\{1+O(\delta(1+x))\},
\qquad 0\le x\le c/\delta,
$$

and $\Lambda(x)=\kappa_3(S)x^3/6+O(\delta^2x^4)$. Therefore the uncorrected normal tail is relatively accurate uniformly for $x=o(\delta^{-1/3})$, and the cubic-corrected tail is relatively accurate for $x=o(\delta^{-1/2})$. Poissonized singleton models give nonunit limiting ratios at each boundary, proving that both exponents are optimal for this class. Both signs and strict/non-strict tails are covered, with endpoint atoms controlled explicitly.

This transfer from cumulants is classical Cramér–Petrov theory; the paper cites Döring–Jansen–Schubert and supplies a self-contained exponential-tilt proof with conservative explicit constants. The contribution here is its application, uniform scope, endpoint control, and sharpness within the batch-discovery model, not a claim to originate the general cumulant method.

### A uniform fixed-sample variance correction

Write $D_{w,n}$ for weighted discovery after exactly $n$ batches, $s_{w,n}^2=\operatorname{Var}(D_{w,n})$, and $f_w(t)=\sum_i w_i\lambda_i e^{-t\lambda_i}$. Uniformly in $n$, the law, and nonnegative weights with batch total at most $W$,

$$
0\le\mathbb ED_{w,n}-\mathbb EX_w(n)\le2W,
\qquad
\left|s_{w,n}^2-\left[V_w(n)-nf_w(n)^2\right]\right|
\le(4+2e^{-2})W^2.
$$

The corrected variance is exactly $\min_a\operatorname{Var}(X_w(n)-aN(n))$. It is asymptotic to the fixed-sample variance whenever either is large compared with $W^2$. The correction orders cannot be improved uniformly over changing laws. A uniform $CW/s_{w,n}$ Kolmogorov bound still needs a distributional proof: the moment correction alone does not establish it.

### Higher coverage with a sharper threshold factor

For integer thresholds $1\le r_i\le r$, let $Y_r=\sum_iw_i\mathbf1_{\{N_i\ge r_i\}}$ and $V_r=\operatorname{Var}(Y_r)$. The exact envelope constant is

$$
\ell_r=1+\sum_{k=1}^{r-1}\frac{e^{-k}k^k}{k!}\sim\sqrt{\frac{2r}{\pi}}.
$$

The higher-coverage paper proves

$$
d_K\!\left(\frac{Y_r-\mathbb EY_r}{\sqrt{V_r}},Z\right)
\le\min\!\left\{1,\frac{15W\sqrt{\ell_r}}{\sqrt{V_r}}\right\}.
$$

This gives $15\sqrt{1+e^{-1}}\,W/\sqrt{V_r}$ at threshold two and $O(r^{1/4}W/\sqrt{V_r})$ for growing thresholds. Disjoint batches with a common threshold give a binomial example attaining the variance order. A star with one high-threshold hub has a fourth cumulant of order $V^{3/2}$ and excludes a threshold-independent variance-linear bound at that order.

The pair paper removes the threshold factor on disjoint unions of stars. More generally, if removing species of total weight at most $H$ in each component leaves independent indicators there, the bound is $C\max(W,H)/\sqrt{V_r}$. On arbitrary pair graphs the proof gives a switch with aggregate quadratic cost at most $3WV_r$ ($2WV_r$ on triangle-free graphs). The open step is to control the signed endpoint bias and the length-weighted crossing probability uniformly in the cutoff. The stars that obstruct the cumulant method satisfy the threshold-uniform Kolmogorov estimate.

## Build and reproduce

The verified Python environment uses Python 3.12 and the package versions pinned in [requirements.txt](requirements.txt). A standard TeX Live installation with `latexmk`, Latin Modern, and the packages imported by the sources builds all six papers; each bibliography is included in its LaTeX file.

```sh
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt
make all
```

For Debian/Ubuntu, the TeX dependencies can be installed with:

```sh
sudo apt-get install latexmk texlive-latex-extra texlive-fonts-recommended lmodern
```

Individual targets:

```sh
make pdf              # Compile all six papers
make check            # Rerun all seven verification suites
make check-optimal    # Check partition, ratio, and cumulant calculations
make check-extensions # Check tails, fixed-sample corrections, and threshold envelopes
make check-fixed      # Check categorical identities, ratios, and sector construction
make check-higher-coverage # Check the threshold-two identity, Stein mass, and capacity calculations
make check-pair-coupling   # Check conditional switches, cumulants, and star forests
make validate         # Build and check labels, citations, logs, PDF text, and margins
make clean            # Remove auxiliary TeX files; preserve PDFs and reports
```

All seven verification scripts accept `--output PATH`. The first five refuse Python optimization mode. The two higher-coverage scripts use explicit checks that remain active in that mode. Each JSON report records the parameters, package versions, and generating script's SHA-256. The checks use deterministic enumeration, finite-state equations, high-precision quadrature, or exact identities; they do not use Monte Carlo estimates of probabilities.

| Foundation checks | Count |
|---|---:|
| Bracket mean and variance cases, thresholds $r=1,2,3$ | 60 |
| Exact threshold means / variance-sandwich cases | 60 each |
| Fixed-sample covariance and correction cases | 25 |
| Centered-indicator moment bounds | 80 |
| Signed-weight variance / fractional MGF / actual-variance MGF | 108 / 144 / 288 |
| Coefficient loads / Hölder parameter-grid cases | 30 / 5,000 |
| Threshold-envelope comparisons / counterexample parameters | 117,160 / 9 |

The optimal-rate suite checks the exact absence/independence-polynomial identity by independent finite enumeration, occupation-ratio bounds on induced subgraphs, weighted and unweighted cumulants, and the numerical smoothing constant. Its [report](reports/optimal_rate_verification.json) records the models, tolerances, case counts, and residuals.

The follow-up suite checks exponential-tilt identities, fixed-sample moments against exact union-state enumeration, the residual-variance projection, threshold envelopes and coefficient loads, proof constants, and binomial boundary examples. Its [report](reports/quantitative_extensions_verification.json) records precision and numerical limitations. These checks supplement the analytic proofs; finite calculations do not establish countable-space or asymptotic theorems.

The fixed-sample suite checks the categorical partition identity against union-state enumeration, induced-subgraph occupation ratios, harmonic parameters, and the two-community coefficient limit. Its [report](reports/fixed_sample_verification.json) distinguishes exact rational identities from finite-precision evaluations; it does not certify the conjectured optimal rate.

The [higher-coverage checks](reports/higher_coverage_verification.json) compare the matching formula with direct threshold enumeration, check two Stein masses, and calculate star cumulants. The [pair-switch checks](reports/pair_coupling_verification.json) compare the zero-bias CDF from switch intervals with the stop-loss formula and account for the omitted Poisson tail.

[Finite-state report](reports/verification_results.json) · [Revision report](reports/revision_verification_results.json) · [Optimal-rate report](reports/optimal_rate_verification.json) · [Follow-up report](reports/quantitative_extensions_verification.json) · [Fixed-sample report](reports/fixed_sample_verification.json)

[Foundation PDF validation](reports/manuscript_validation.json) · [Optimal-rate PDF validation](reports/optimal_manuscript_validation.json) · [Follow-up PDF validation](reports/quantitative_extensions_validation.json) · [Fixed-sample PDF validation](reports/fixed_sample_validation.json) · [Higher-coverage PDF validation](reports/higher_coverage_validation.json) · [Pair-switch PDF validation](reports/pair_indicator_coupling_validation.json)

The [GitHub Actions workflow](.github/workflows/verify.yml) builds all six papers, runs all seven suites and the document checks, and uploads the PDFs and reports as workflow artifacts. Visual page review remains a separate manual check.

## Provenance and scope

The complete November 12 revised bundle is preserved byte-for-byte in [the archive](archive/2026-11-12/README.md), including its original source, PDF, scripts, audit notes, and reports. The foundation source changes are documented in the [initial November 26 audit](audit/2026-11-26.md). The companion theorem and subsequent results are recorded separately in the [optimal-rate audit](audit/optimal_rate_2026-11-26.md), [quantitative-extensions audit](audit/quantitative_extensions_2026-11-26.md), [fixed-sample audit](audit/fixed_sample_2026-11-26.md), and [higher-coverage audit](audit/higher_coverage_and_pairs_2026-11-27.md). Earlier manuscripts and audit records retain their historical scope statements; the current status is given above and in the [research-status note](extensions/research_status.md).

The optimal rate and relative tails apply to scalar Poissonized first-discovery counts and their nonnegative weighted versions. The fixed-sample paper supplies a rate with a logarithmic factor, also for signed weights, and leaves its removal open. A variance-only signed-weight CLT, uniform fixed-sample moderate deviations and relative tails, joint occupancy-spectrum and functional limits, local limits with appropriate arithmetic hypotheses, and unknown-law prediction intervals require further work. The threshold-uniform higher-coverage rate for unrestricted pair graphs remains open. A finite second moment alone remains insufficient for the unbounded-batch CLT, as shown by the foundation counterexample.

## Citation and licensing

Use each manuscript's author, title, and revision date when citing these working drafts; [BibTeX entries](CITATION.bib) are provided. No journal publication or DOI is claimed for these drafts.

The [MIT License](LICENSE) applies to the verification scripts and repository build/CI infrastructure. The mathematical manuscripts, compiled PDFs, research notes, and archived research materials are excluded from that grant; this setup makes no change to their licensing. Authorship of the manuscripts remains Christopher D. Long.
