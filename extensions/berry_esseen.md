# Berry–Esseen bounds beyond sparse dependency graphs

Research note — November 26, 2026

## Background and status

The manuscript [A fluctuation theory for bounded-batch species discovery](../papers/bounded_batch_discovery.pdf), first drafted on November 12 and audited for this repository on November 26, proves a variance-normalized CLT with a uniform quantitative rate for Poissonized bounded-batch discovery, without a bound on the maximum degree of the species dependency graph.

Maximilian Janisch's [A truncated-moment Berry–Esseen bound for dependency graphs](https://github.com/maximilian-janisch/berry-esseen-dependency-graphs/blob/8a2c9e9a4eee068413f0a5fdfbf8e16172870c8a/paper/janisch-lehericy-proof.pdf) reports a proof of the Raič and Janisch–Lehéricy conjectures. This note uses Theorem 1 and its centered-moment Corollary 2 at commit [`8a2c9e9`](https://github.com/maximilian-janisch/berry-esseen-dependency-graphs/commit/8a2c9e9a4eee068413f0a5fdfbf8e16172870c8a), dated November 26, 2026. The theorem statements and relevant formal definitions were inspected; the Lean formalization has not been independently built or replayed for this project. The application below is conditional on that external theorem. See the [reference audit](../audit/references_2026-11-26.md).

## A bounded-degree corollary

Let $A$ be a random finite batch, $|A|\le B$, arriving according to a rate-one Poisson process. Put
$$
\lambda_i=\Pr(i\in A),\quad
\lambda_{ij}=\Pr(i,j\in A),\quad
I_i(t)=\mathbf{1}_{\{i\text{ discovered by }t\}},\quad
p_i=1-e^{-t\lambda_i}.
$$
Write $X=\sum_i I_i(t)$ and $V=\operatorname{Var}(X)>0$.

Join distinct species $i,j$ if $\lambda_{ij}>0$. Poisson splitting into independent batch-type processes proves that the resulting graph is a dependency graph for the indicators: if no edge joins two sets of species, no positive-rate batch meets both sets, so their indicator vectors depend on disjoint independent coordinate families. If the maximum degree is $D<\infty$, Janisch's centered truncated-moment bound, in its third-moment specialization, gives
$$
d_K((X-\mathbb EX)/\sqrt V,Z)
\le \min\left\{1,\frac{250000(D+1)^2}{V^{3/2}}\sum_i\mathbb E|I_i-p_i|^3\right\},
$$
where $Z$ is standard normal.

Since
$$
\mathbb E|I_i-p_i|^3=p_i(1-p_i)[p_i^2+(1-p_i)^2]\le p_i(1-p_i),
$$
and discovery indicators have nonnegative covariances (by their independent-clock representation), the sum of third moments is at most $V$. Consequently,
$$
d_K((X-\mathbb EX)/\sqrt V,Z)
\le \min\{1,250000(D+1)^2/\sqrt V\}.
$$
For countably many species, apply the finite result to increasing finite restrictions. The counts converge in $L^2$ by domination by $BN(t)$, so means and variances converge and the continuous normal comparison passes to the limit.

## Why bounded batches are different

The overlap identity
$$
\sum_{j\ne i}\lambda_{ij}\le (B-1)\lambda_i
$$
holds regardless of the graph's degree. Uniformly sampled pairs on $n$ species have $B=2$, $D=n-1$, and $\lambda_{ij}/\lambda_i=1/(n-1)$ for distinct species. Thus a bound involving $D$ can be vacuous even in a bounded-batch model.

A weighted dependency theorem must account for the strength and joint structure of co-discovery; a formal replacement of degree by a sum of pairwise weights in an arbitrary dependency-graph theorem is not justified.

## Open proof target

**Conjecture.** There exists a universal constant $C$ such that, for all Poissonized bounded-batch laws and times with $V>0$,
$$
d_K((X-\mathbb EX)/\sqrt V,Z)
\le \min\{1,CB/\sqrt V\}.
$$

This is an open target here, not an asserted new theorem or an exhaustive novelty claim. For a concrete obstruction, choose each of $N$ disjoint packages of $B$ species with probability $1/N$, and observe at time $t=cN$ for fixed $c>0$. Then $X=B\,\mathrm{Bin}(N,p)$ with $p=1-e^{-c}$ and $V=B^2Np(1-p)$. The largest atom has order $N^{-1/2}$, equivalently $B/\sqrt V$, by the binomial local limit estimate. Since any continuous comparison distribution differs from a distribution with an atom of size $a$ by Kolmogorov distance at least $a/2$, this shows the proposed order cannot be improved uniformly.

Research directions:

- Adapt the Stein argument to the independent batch coordinates and the weighted overlap identity.
- Combine bounded changes with the archived manuscript's predictable-bracket variance estimate.
- Investigate size-bias and Poisson difference-operator couplings, and determine whether existing normal approximation theorems already give the conjectured rate.
- Extend any optimal-order result to nonnegative weighted discovery and multiple-coverage counts, then investigate quantitative de-Poissonization.

This note is separate from the manuscript and supplies no new proof of the conjecture. The complete original November 12 bundle is preserved unchanged in [the archive](../archive/2026-11-12/README.md); the current manuscript contains only the corrections and reproducibility updates recorded in the [repository audit](../audit/2026-11-26.md).
