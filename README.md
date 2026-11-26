# Bounded-batch species discovery

**A fluctuation theory for bounded-batch species discovery**

Christopher D. Long — research draft, November 12, 2026; repository revision, November 26, 2026.

[Read the paper (PDF)](papers/bounded_batch_discovery.pdf) · [LaTeX source](papers/bounded_batch_discovery.tex) · [Audit](audit/2026-11-26.md) · [Original bundle](archive/2026-11-12/README.md)

This project studies the number of distinct species discovered through independent random batches. Species within a batch may be dependent. A batch may be a graph edge, a bounded hyperedge, or another finite set of observations.

The current manuscript contains proofs of the results below. A fresh internal mathematical and reference audit found no substantive error in the arguments reviewed and made small domain and degenerate-case clarifications. This is a research draft, not independent peer review or formal verification. The [audit record](audit/2026-11-26.md) describes the review and its limits.

## Main results

For i.i.d. batches $A_1,A_2,\ldots$ on a countable species space, with a **fixed deterministic law** and $|A_k|\le B$, write

$$
D_n=\left|\bigcup_{k=1}^{n}A_k\right|,\qquad
X(t)=D_{N(t)},\qquad V(t)=\operatorname{Var}X(t),
$$

where $N$ is an independent rate-one Poisson process. Set $\lambda_i=\Pr(i\in A_1)$ and $v_0(t)=\sum_i e^{-t\lambda_i}(1-e^{-t\lambda_i})$.

| Result | Scope and manuscript location |
|---|---|
| Variance-normalized central limit theorem whenever $\operatorname{Var}(D_n)\to\infty$ | Fixed law; no connectivity, regular variation, or within-batch independence assumption. Theorem 1.1 and Section 7. |
| Sharp variance comparison $v_0(t)\le V(t)\le Bv_0(t)$ | Poissonized discoveries. Section 3. |
| Uniform quantitative normal approximation | Poissonized bounded-batch laws, including changing laws. Section 5. |
| Actual-variance Bernstein concentration and Gaussian moderate deviations | Poissonized counts. Section 9. |
| Absolute expansion $\operatorname{Var}(D_n)=V(n)-n m'(n)^2+o(1)$ | Fixed laws with $\mathbb E|A_1|^2<\infty$, where $m(t)=\mathbb E X(t)$. Section 6. |
| Lindeberg extension and a finite-second-moment counterexample | Unbounded batches; finite second moment alone does not imply a CLT. Section 8. |
| Threshold, nonnegative weighted, conditional, and temporal results | Each with its stated assumptions. Sections 10–11. |

With $\eta=B^2/V(t)$, the Poissonized Kolmogorov bound is

$$
d_K\!\left(\frac{X(t)-\mathbb E X(t)}{\sqrt{V(t)}},Z\right)
\le \min\!\left\{1,C\left[\eta\log(e+\eta^{-1})\right]^{1/4}\right\},
\qquad Z\sim\mathcal N(0,1).
$$

The manuscript also gives the simpler $C\eta^{1/5}$ bound and

$$
\Pr\bigl(|X(t)-\mathbb E X(t)|\ge x\bigr)
\le 2\exp\!\left(-\frac{x^2}{2(V(t)+2Bx)}\right),\qquad x>0.
$$

For $B=2$, the fixed-law theorem addresses the vertex-count question in Janson's countable-edge model without the connectivity hypothesis used in Eriksson's connected-regime result. The [reference audit](audit/references_2026-11-26.md) checks the cited problem and theorem statements. Diffuse dust mechanisms and arbitrary mixtures over random laws are outside this claim.

## Berry–Esseen research note

[Berry–Esseen bounds beyond sparse dependency graphs](extensions/berry_esseen.md) applies the centered-moment corollary stated in Maximilian Janisch's November 26, 2026 manuscript to discovery models with maximum dependency degree $D<\infty$. It records the conditional bound $\min\{1,250000(D+1)^2/\sqrt V\}$, pins the inspected upstream commit, and explains why bounded batch size does not bound dependency degree.

The general bound $CB/\sqrt V$ remains an **open proof target in this repository**. The note is separate from the manuscript's proofs. Janisch's Lean formalization was inspected at its theorem boundary but was not independently replayed here.

## Build and reproduce

The verified Python environment uses Python 3.12 and the package versions pinned in [requirements.txt](requirements.txt). A standard TeX Live installation with `latexmk`, Latin Modern, and the packages imported by the source builds the paper; the bibliography is included in the LaTeX file.

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
make pdf       # Compile papers/bounded_batch_discovery.pdf
make check     # Rerun both numerical/symbolic verification suites
make validate  # Build and check labels, citations, log, PDF text, and margins
make clean     # Remove auxiliary TeX files; preserve the PDF and reports
```

Both verification scripts also accept `--output PATH`. They refuse Python optimization mode, which would disable their assertions. Each JSON report records the parameters, package versions, and generating script's SHA-256. The finite-state suite uses a fixed seed to select two example laws; computations enumerate states or solve finite-state equations, with no Monte Carlo sampling.

| Current executed checks | Count |
|---|---:|
| Bracket mean and variance cases, thresholds $r=1,2,3$ | 60 |
| Exact threshold means / variance-sandwich cases | 60 each |
| Fixed-sample covariance and correction cases | 25 |
| Centered-indicator moment bounds | 80 |
| Signed-weight variance / fractional MGF / actual-variance MGF | 108 / 144 / 288 |
| Coefficient loads / Hölder parameter-grid cases | 30 / 5,000 |
| Threshold-envelope comparisons / counterexample parameters | 117,160 / 9 |

All passed in the November 26 audit. The maximum absolute discrepancy in $\mathbb E Q=V$ was $6.26\times10^{-11}$. The added $r=3$ cases exercise the improved envelope $d_3=3e^2<32$. These finite numerical and symbolic checks supplement the proofs; they do not establish asymptotic or countable-space theorems.

[Finite-state report](reports/verification_results.json) · [Revision report](reports/revision_verification_results.json) · [Manuscript validation](reports/manuscript_validation.json)

The [GitHub Actions workflow](.github/workflows/verify.yml) builds the paper, runs both suites and the document checks, and uploads the PDF and reports as workflow artifacts. Visual page review remains a separate manual check.

## Provenance and scope

The complete November 12 revised bundle is preserved byte-for-byte in [the archive](archive/2026-11-12/README.md), including its original source, PDF, scripts, audit notes, and reports. The current source changes are documented in the [November 26 audit](audit/2026-11-26.md). Historical validation records retain their original dates.

No optimal $O(B/\sqrt V)$ general rate, variance-only fixed-sample rate, fixed-sample moderate-deviation principle, relative Gaussian-tail approximation, joint occupancy-spectrum CLT, functional CLT, or unknown-law prediction interval is claimed. Weighted discoveries use the bracket-load bound $2W^2$ for a positive upper bound $W$ on total batch weight.

## Citation and licensing

Use the manuscript's author, title, and revision date when citing this working draft; a [BibTeX entry](CITATION.bib) is provided. No journal publication or DOI is claimed.

The [MIT License](LICENSE) applies to the verification scripts and repository build/CI infrastructure. The mathematical manuscript, compiled PDF, research notes, and archived research materials are excluded from that grant; this setup makes no change to their licensing. Authorship of the manuscript remains Christopher D. Long.
