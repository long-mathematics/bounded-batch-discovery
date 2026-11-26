# Bounded-batch species discovery

Research by Christopher D. Long on fluctuations in dependent species-discovery and occupancy models.

## Manuscript

[A fluctuation theory for bounded-batch species discovery](papers/bounded_batch_discovery.tex), dated November 12, 2026 (revised research draft; proof and reference audit).

The original LaTeX source is archived separately from later research extensions. The draft establishes, subject to independent mathematical review:

- A central limit theorem for the number of discovered species under a fixed law of independent bounded batches, whenever the actual variance diverges.
- The countable-edge vertex-count case of Janson's question, without the connectivity hypothesis.
- Sharp variance comparisons, quantitative Poissonized normal approximation, concentration, and moderate deviations.
- A fixed-sample de-Poissonization variance expansion, a Lindeberg extension, and results for multiple-coverage and weighted counts.

For Poissonized variance \(V\) and maximum batch size \(B\), the manuscript gives a Kolmogorov bound of order
\[
\left[\frac{B^2}{V}\log\left(e+\frac{V}{B^2}\right)\right]^{1/4}.
\]

## Current research direction

[Berry–Esseen bounds beyond sparse dependency graphs](extensions/berry_esseen.md) records the connection to Maximilian Janisch's November 26, 2026 manuscript and the proposed optimal-order bound for bounded-batch discovery, including models with dense dependency graphs.

The November 12 manuscript is a research draft, not a claim of independent peer review or formal verification. The proposed optimal-order bound is an open proof target in this repository, not an established theorem.

## Build

With a standard LaTeX distribution containing the packages listed in the source:

```sh
cd papers
latexmk -pdf bounded_batch_discovery.tex
```

The original PDF is not yet included in this repository.
