# Berry–Esseen bounds beyond sparse dependency graphs

Research note — November 26, 2026; updated after the optimal-rate proof and quantitative follow-up.

## Background and status

**The general target is now proved in the companion paper [Optimal Berry–Esseen bounds for bounded-batch discovery](../papers/optimal_berry_esseen.pdf):**

$$
d_K((X-\mathbb EX)/\sqrt V,Z)\le\min\{1,150B/\sqrt V\}.
$$

The proof allows arbitrary overlap and countably many species, and extends to nonnegative weights with total contribution at most $W$ per batch. It uses a bipartite independence-polynomial representation, a zero-free disk, the foundation paper's actual-variance moment bound, and Fourier smoothing. It does not rely on Janisch's dependency-graph theorem. See the [source](../papers/optimal_berry_esseen.tex) and [proof audit](../audit/optimal_rate_2026-11-26.md).

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

## Resolution of the original proof target

**Theorem (companion paper, Theorems 1.1–1.2).** For all Poissonized bounded-batch laws and times with $V>0$,
$$
d_K((X-\mathbb EX)/\sqrt V,Z)
\le \min\{1,150B/\sqrt V\}.
$$

The same statement holds with $B$ replaced by $W$ for weighted discoveries, where $W>0$ bounds the total nonnegative weight of each batch; a cardinality bound is unnecessary. The proof gives the stronger numerical coefficient $36/\sqrt\pi+864/(\pi\sqrt{2\pi})<131$, rounded upward to 150 in the main theorem.

For the key zero-free step, form the bipartite incidence graph between species and positive-probability batch types. Give species $i$ activity $u_i$ and batch type $C$ activity $e^{t\mu(C)}-1$. The absence polynomial is a positive constant times this graph's independence polynomial. An induction on induced subgraphs proves nonvanishing whenever

$$
|u_i|\le r_i<1,\qquad
\sum_{i\in C}\arcsin r_i<\pi/2.
$$

The species occupation ratios remain inside their disks, while the batch ratios remain in the closed right half-plane. Species degrees never enter the estimate. Substituting $u_i=e^{-zh_i}-1$ gives a zero-free weighted MGF for $|z|<1/(2W)$.

Combine this with the already proved real-axis estimate

$$
\log\mathbb E e^{\theta(X_h-\mathbb EX_h)}
\le \frac{V_h\theta^2}{2(1-2W|\theta|)}.
$$

Harmonic coefficient estimates on $|z|=1/(4W)$ give $|\kappa_k(X_h)|\le2k!V_h(4W)^{k-2}$. Summing the cumulant remainder and applying Esseen's smoothing inequality proves the optimal rate. The variance-linear estimate uses both zero-freeness and the real-axis moment bound; zero-freeness alone would not suffice.

The order cannot be improved uniformly. Choose each of $N$ disjoint packages of $B$ species with probability $1/N$, and observe at time $t=N\log2$. Then $X=B\operatorname{Bin}(N,1/2)$ and $V=B^2N/4$. Its largest atom has order $N^{-1/2}$, equivalently $B/\sqrt V$. Any continuous comparison distribution has Kolmogorov distance at least half the largest atom. The companion paper states this lower bound explicitly.

## Remaining questions

The follow-up [Further quantitative bounds for bounded-batch discovery](../papers/quantitative_extensions.pdf) now proves relative Gaussian tails, a cubic correction, and sharpness of their ranges. It also proves a uniform $O(W^2)$ fixed-sample variance correction and improves the higher-coverage bracket factor from order $r$ to order $\sqrt r$.

The optimal fixed-sample and higher-coverage rates remain open here. Quantitative de-Poissonization needs a distributional argument at the corrected variance scale. Higher coverage needs an analytic or probabilistic replacement for the first-discovery factorization; the straightforward truncated-Poisson factorization is false. The ordinary Poisson add-one Stein mass also has an unavoidable threshold-dependent variance in a singleton example, which obstructs that proposed intermediate estimate without disproving the conjectured rate. The [research-status note](research_status.md) records these precise gaps and examples. None of the rejected Stein estimates is used in the completed first-discovery proof.

The complete original November 12 bundle is preserved unchanged in [the archive](../archive/2026-11-12/README.md). The foundation manuscript retains the corrections recorded in the [initial repository audit](../audit/2026-11-26.md); the new theorem is in a separate companion manuscript. No exhaustive novelty claim is made for the auxiliary zero-free lemma or the analytic tools.
