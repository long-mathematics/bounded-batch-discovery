# Bounded-batch species discovery

**Optimal Berry–Esseen bounds, sharp relative-tail windows, and quantitative extensions for discovery through random batches.**

Christopher D. Long — research drafts, November 2026.

[Further quantitative bounds (PDF)](papers/quantitative_extensions.pdf) · [LaTeX source](papers/quantitative_extensions.tex) · [Proof audit](audit/quantitative_extensions_2026-11-26.md) · [Research status and open questions](extensions/research_status.md)

[Optimal-rate paper (PDF)](papers/optimal_berry_esseen.pdf) · [LaTeX source](papers/optimal_berry_esseen.tex) · [Proof audit](audit/optimal_rate_2026-11-26.md)

[Foundation paper (PDF)](papers/bounded_batch_discovery.pdf) · [LaTeX source](papers/bounded_batch_discovery.tex) · [Foundation audit](audit/2026-11-26.md) · [Original bundle](archive/2026-11-12/README.md)

This project studies the number of distinct species discovered through independent random batches. Species within a batch may be dependent. A batch may be a graph edge, a bounded hyperedge, or another finite set of observations.

The November 26 optimal-rate paper proves the original target in this repository: a uniform Poissonized Kolmogorov bound of order $B/\sqrt V$, with the explicit constant 150. The follow-up paper proves relative Gaussian tails and a cubic correction with sharp ranges, uniform fixed-sample mean and variance corrections, and sharper bounds for higher coverage. The optimal fixed-sample and higher-coverage Berry–Esseen questions remain open here. Internal mathematical audits checked the arguments and their imported estimates; the papers remain research drafts, without independent peer review or formal verification. The linked audit records describe the checks and their limits.

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
| Threshold bracket load $2\ell_rW^2$, where $\ell_r\sim\sqrt{2r/\pi}$ | Nonnegative weights and heterogeneous integer thresholds at most $r$; improves generic normal approximation and actual-variance concentration. Follow-up Section 8. |
| Variance-normalized central limit theorem whenever $\operatorname{Var}(D_n)\to\infty$ | Fixed law; no connectivity, regular variation, or within-batch independence assumption. Foundation Theorem 1.1 and Section 7. |
| Sharp variance comparison $v_0(t)\le V(t)\le Bv_0(t)$ | Poissonized discoveries. Foundation Section 3. |
| Actual-variance Bernstein concentration and Gaussian moderate deviations | Poissonized counts. Foundation Section 9. |
| Absolute expansion $\operatorname{Var}(D_n)=V(n)-n m'(n)^2+o(1)$ | Fixed laws with $\mathbb E|A_1|^2<\infty$, where $m(t)=\mathbb E X(t)$. Foundation Section 6. |
| Lindeberg extension and a finite-second-moment counterexample | Unbounded batches; finite second moment alone does not imply a CLT. Foundation Section 8. |
| Threshold, conditional, and temporal results | Each with its stated assumptions. Foundation Sections 10–11; the follow-up improves threshold bounds. The optimal rate concerns first discovery. |

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

The generic martingale normal-approximation bound improves to

$$
d_K\!\left(\frac{Y_r-\mathbb EY_r}{\sqrt{V_r}},Z\right)
\le\min\!\left\{1,C\left[\eta_r\log(e+\eta_r^{-1})\right]^{1/4}\right\},
\qquad\eta_r=\frac{\ell_rW^2}{V_r}.
$$

Thus $\sqrt r\,W^2/V_r\to0$ suffices for a CLT, including changing laws and growing thresholds. Actual-variance Bernstein concentration has linear scale $2W\sqrt{\ell_r}=O(Wr^{1/4})$. The constant $\ell_r$ is sharp for the deterministic envelope inequality; this does not establish that threshold dependence is necessary in a probabilistic bound. The conjectural threshold-uniform $CW/\sqrt{V_r}$ rate remains open here.

## Build and reproduce

The verified Python environment uses Python 3.12 and the package versions pinned in [requirements.txt](requirements.txt). A standard TeX Live installation with `latexmk`, Latin Modern, and the packages imported by the sources builds all three papers; each bibliography is included in its LaTeX file.

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
make pdf              # Compile all three papers
make check            # Rerun all four verification suites
make check-optimal    # Check partition, ratio, and cumulant calculations
make check-extensions # Check tails, fixed-sample corrections, and threshold envelopes
make validate         # Build and check labels, citations, logs, PDF text, and margins
make clean            # Remove auxiliary TeX files; preserve PDFs and reports
```

All four verification scripts accept `--output PATH`. They refuse Python optimization mode, which would disable their assertions. Each JSON report records the parameters, package versions, and generating script's SHA-256. The checks use deterministic enumeration, finite-state equations, high-precision quadrature, or exact identities; they do not use Monte Carlo estimates of probabilities.

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

[Finite-state report](reports/verification_results.json) · [Revision report](reports/revision_verification_results.json) · [Optimal-rate report](reports/optimal_rate_verification.json) · [Follow-up report](reports/quantitative_extensions_verification.json)

[Foundation PDF validation](reports/manuscript_validation.json) · [Optimal-rate PDF validation](reports/optimal_manuscript_validation.json) · [Follow-up PDF validation](reports/quantitative_extensions_validation.json)

The [GitHub Actions workflow](.github/workflows/verify.yml) builds all three papers, runs all four suites and the document checks, and uploads the PDFs and reports as workflow artifacts. Visual page review remains a separate manual check.

## Provenance and scope

The complete November 12 revised bundle is preserved byte-for-byte in [the archive](archive/2026-11-12/README.md), including its original source, PDF, scripts, audit notes, and reports. The foundation source changes are documented in the [initial November 26 audit](audit/2026-11-26.md). The companion theorem and follow-up results are recorded separately in the [optimal-rate audit](audit/optimal_rate_2026-11-26.md) and [quantitative-extensions audit](audit/quantitative_extensions_2026-11-26.md). Earlier manuscripts and audit records retain their historical scope statements; the current status is given above and in the [research-status note](extensions/research_status.md).

The optimal rate and new relative tails apply to scalar Poissonized first-discovery counts and their nonnegative weighted versions. The new fixed-sample moment corrections and higher-coverage bounds do not resolve optimal rates for those models. Signed weights, uniform fixed-sample moderate deviations and relative tails, joint occupancy-spectrum and functional limits, local limits with appropriate arithmetic hypotheses, and unknown-law prediction intervals require further work. A finite second moment alone remains insufficient for the unbounded-batch CLT, as shown by the foundation counterexample.

## Citation and licensing

Use each manuscript's author, title, and revision date when citing these working drafts; [BibTeX entries](CITATION.bib) are provided. No journal publication or DOI is claimed for these drafts.

The [MIT License](LICENSE) applies to the verification scripts and repository build/CI infrastructure. The mathematical manuscripts, compiled PDFs, research notes, and archived research materials are excluded from that grant; this setup makes no change to their licensing. Authorship of the manuscripts remains Christopher D. Long.
