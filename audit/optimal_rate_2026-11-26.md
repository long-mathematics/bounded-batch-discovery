# Proof audit of the optimal Poissonized Berry–Esseen rate

Date: November 26, 2026.

## Assessment and reviewed sources

The argument in [Optimal Berry–Esseen bounds for bounded-batch discovery](../papers/optimal_berry_esseen.tex) supports the stated uniform bound

$$
d_{\mathrm K}\!\left(\frac{X-\mathbb EX}{\sqrt V},\mathcal N(0,1)\right)
\le \min\left\{1,\frac{150B}{\sqrt V}\right\}.
$$

The same proof supports its nonnegative weighted version with total batch weight at most $W$, replacing $B$ by $W$. No unresolved mathematical gap was identified in the proof steps reviewed below. This is a research proof audit; it is not formal verification, independent peer review, or a determination of novelty or priority.

The source examined for this audit has SHA-256
`0084f5b37d2be440d2e5762520055c5c9b43787aa342dfa6e50c6307d4ab860c`.
The imported moment-generating-function argument was rechecked in the
[companion manuscript](../papers/bounded_batch_discovery.tex), whose source at the time of this check has SHA-256
`360523998d267ca8226d41dfeaeeb33b1a32b8b414cdc5c867034bdd45f4a84b`.
These fingerprints identify the final reviewed versions, including the explicit restriction to positive-probability batch types and the integer cardinality in the sharpness example. Later editorial changes may have different hashes.

The proof combines two distinct ingredients: a complex zero-free disk proved in the new note, and a real-axis exponential-moment bound with the actual variance, imported from the companion manuscript. Neither ingredient is silently substituted for the other.

## Model and exact scope

The species space is countable. The batches are independent and identically distributed random **finite** subsets of that space, with a deterministic sampling law. They arrive according to an independent rate-one Poisson process. Empty batches are allowed. At a fixed finite time $t$, $X(t)$ counts species seen at least once.

For the unweighted theorem, $|A_1|\le B$ almost surely, with $B\in\mathbb N$. For the weighted theorem, deterministic weights $h_i\ge0$ satisfy

$$
\sum_{i\in A_1}h_i\le W\quad\text{almost surely},\qquad 0<W<\infty.
$$

The weighted theorem does not require a uniform bound on batch cardinality, although each individual batch remains finite. Species of weight zero or zero marginal arrival rate may be removed. Every remaining species then has $0<h_i\le W$. No assumption is made on a species' number of neighbors, the density of overlap, connectedness, or independence within a batch.

The normalized conclusions require $V>0$. If the variance is zero, the centered count is zero almost surely and the displayed normalization is undefined; the theorem correctly excludes that case. The domination

$$
0\le X_h(t)\le WN(t)
$$

gives every exponential moment used in the argument. The bound is uniform in the law and all other parameters, so it applies to sequences in which the law, weights, time, and species space vary.

The new theorem concerns scalar **Poissonized first-discovery counts**. It does not prove the same optimal rate for a deterministic sample size, counts seen at least $r\ge2$ times, signed weights, joint distributions, a local limit theorem, or relative Gaussian tails. The note states these limits of its scope explicitly.

## 1. Bipartite partition-function lemma

The finite-graph induction was checked for every induced subgraph, including isolated vertices and zero batch activities. The left vertices carry complex activities $z_i$ with $|z_i|\le r_i<1$; the right vertices carry nonnegative real activities. The condition is

$$
\sum_{i\in N(j)}\arcsin r_i<\frac\pi2
$$

for each right vertex $j$.

For a vertex $v$, the occupation ratio is

$$
R_v(H)=\lambda_v\frac{\mathcal Z_{H-v-N_H(v)}}{\mathcal Z_{H-v}}.
$$

The exact product recursion deletes $v$, then deletes its neighbors successively. Consequently every ratio used on its right side belongs to a **proper** induced subgraph. The induction therefore justifies all denominators before using them; it is not circular.

For a left vertex, every neighboring right ratio has nonnegative real part, so the factors $|1+R_j|$ are at least one. This proves $|R_i|\le r_i$. For a right vertex, the elementary disk geometry

$$
|\arg(1+\zeta)|\le\arcsin r\qquad(|\zeta|\le r<1)
$$

keeps the argument of the product of inverse factors strictly between $-\pi/2$ and $\pi/2$. Multiplication by a positive real activity preserves this interval, giving a positive real part; zero activity gives ratio zero. Finally $1+R_v\ne0$ on either side, so the deletion identity proves nonvanishing of the partition function.

The one-sided degree corollary follows with radius $\sin(\pi/(2B))$. It imposes no bound on degrees on the other side. For $B=1$, the strict activity inequality permits a common auxiliary radius below one, as required by the lemma.

## 2. Finite discovery and the zero-free MGF

Restriction to a finite species set is legitimate even when the original batch law has countably many supported types. Intersect the batches with the finite set, discard empty intersections, and coalesce equal resulting subsets. Independent Poisson thinning then supplies finitely many independent batch-type counts. Their total intensity is finite, and the cardinality or weight bound is preserved.

Writing $a_C$ for a resulting type's Poisson mean and $q_C=e^{-a_C}$, the absence polynomial is exactly

$$
Z_J(\boldsymbol u)
=\sum_{S\subseteq J}\left(\prod_{i\in S}u_i\right)
  \prod_{C:C\cap S\ne\varnothing}q_C.
$$

Assigning batch activity $q_C^{-1}-1\ge0$ in the species–batch incidence graph gives

$$
Z_J(\boldsymbol u)=\left(\prod_Cq_C\right)\mathcal Z_G.
$$

This normalization was checked by summing the permitted batch vertices after fixing the chosen species vertices. The prefactor is positive and cannot introduce zeros.

For weighted discovery,

$$
\mathbb E e^{zX_{h,J}}
=e^{z\sum_{i\in J}h_i}
 Z_J\bigl((e^{-zh_i}-1)_{i\in J}\bigr).
$$

For $|z|<1/(2W)$, choose $r_i=e^{|z|h_i}-1$. Since $|z|h_i<1/2$, the inequality $e^x-1\le2x$ yields $r_i<1$ and

$$
\sum_{i\in C}r_i\le2|z|\sum_{i\in C}h_i<1.
$$

Convexity of $\arcsin$ on $[0,1]$ gives its chord bound $\arcsin r\le(\pi/2)r$. Thus the angular hypothesis holds for every restricted type. This proves the claimed zero-free MGF disk $|z|<1/(2W)$, independently of batch cardinality. The separate unweighted radius $\log(1+\sin(\pi/(2B)))$ also follows from the same representation.

## 3. Countably many species

For an increasing finite exhaustion, put $F_n=X_{h,J_n}$ and $F=X_h(t)$. Then $0\le F_n\uparrow F\le WN(t)$. The estimate

$$
\sup_{|z|\le A}|e^{zF_n}-e^{zF}|
\le A(F-F_n)e^{AF}
$$

follows by integrating the derivative of $x\mapsto e^{zx}$ between $F_n$ and $F$. Its upper bound $AF e^{AF}$ is integrable because $N(t)$ is Poisson. Dominated convergence therefore proves locally uniform convergence of the MGFs in the entire complex plane, not merely pointwise convergence on the real axis.

The limit is entire. Hurwitz's theorem on the connected disk $|z|<1/(2W)$ gives a zero-free limit because its value at zero is one. Centering multiplies it by a nonvanishing exponential. Means and variances also converge by the same Poisson domination. No finite total sum of all species weights is assumed or needed.

An alternative countable passage would apply the finite-model Kolmogorov bound and use weak convergence plus continuity of the normal distribution function. The note's stronger MGF passage directly supports both its cumulant and distributional conclusions.

## 4. Imported actual-variance MGF bound

The dependency on Proposition 9.1 and Theorem 11.1 of the companion manuscript was rechecked, including their underlying bracket calculation. The relevant terminal martingale has jumps

$$
H_C(s)=\sum_{i\in C}h_i e^{-\lambda_i(t-s)}
       \mathbf1_{\{i\text{ absent just before }s\}},
\qquad 0\le H_C(s)\le W.
$$

Its predictable quadratic variation $Q_t$ has mean $V$ and a nonnegative decomposition into functions between zero and one. The weighted incidence identity

$$
\sum_jh_j\lambda_{ij}
=\mathbb E\!\left[\mathbf1_{\{i\in A_1\}}
                  \sum_{j\in A_1}h_j\right]
\le W\lambda_i
$$

bounds each independent batch-coordinate load by $2W^2$. The companion product-space exponential inequality consequently gives

$$
\log\mathbb E e^{a(Q_t-V)}
\le\frac{V}{2W^2}\bigl(e^{2W^2a}-1-2W^2a\bigr).
$$

The compensated-Poisson exponential martingale is a true martingale: the integrands are bounded, the total event rate is finite, and finite-horizon stopped versions have an integrable Poisson exponential bound. Its compensator is at most $\psi_W(|\theta|)Q_t$, where

$$
\psi_W(u)=\frac{e^{Wu}-1-Wu}{W^2}.
$$

Hölder's inequality with exponents $q$ and $q/(q-1)$ therefore produces the displayed bound in the new note. It is not an assumption of independence between the terminal martingale and its bracket.

With $b=W|\theta|\in(0,1/2)$ and $q=1+b$, the two series estimates yield

$$
\log\mathbb E e^{\theta M}
\le\frac{V\theta^2}{2}\frac{q}{1-qb/3-q^2b/2}.
$$

The denominator is positive, and comparison to $q(1-2b)$ follows from the exact identity

$$
1-\frac{qb}{3}-\frac{q^2b}{2}-q(1-2b)
=\frac b6+\frac{2b^2}{3}-\frac{b^3}{2}\ge0.
$$

This proves the imported bound with the actual variance:

$$
\log\mathbb E e^{\theta M}
\le\frac{V\theta^2}{2(1-2W|\theta|)}.
$$

The positive bracket decomposition and its weight calculation remain valid under the weighted hypothesis without a cardinality bound. Countable limits are justified by the bounded-jump and Poisson domination arguments already present in the companion proof.

## 5. Why zero-freeness alone would not suffice

A fixed-size zero-free MGF disk does not imply a universal bound on higher cumulants in terms of the variance alone. For example, if $F_n=n\,\mathrm{Bernoulli}(e^{-n})$, its MGF zeros satisfy

$$
z=\frac{\log((1-e^{-n})/e^{-n})+(2k+1)\pi i}{n}.
$$

Their distances from zero stay bounded away from zero, while
$\kappa_3(F_n)/\operatorname{Var}(F_n)=n(1-2e^{-n})$ diverges. This is not a counterexample to the discovery result with a fixed batch bound. It identifies why the actual-variance real-axis estimate is an essential proof input. The new note correctly makes that dependency explicit.

## 6. Harmonic coefficient estimate and smoothing constant

The centered MGF has an analytic logarithm $K$ on the simply connected zero-free disk, normalized by $K(0)=0$. Let $R=1/(4W)$. The real-axis estimate and the modulus inequality give

$$
\Re K(z)\le VR^2\qquad(|z|\le R).
$$

There is no boundary regularity issue: this closed circle lies strictly inside the larger analytic disk of radius $1/(2W)$.

For $K(z)=\sum_{k\ge1}c_kz^k$, the nonnegative function
$VR^2-\Re K(Re^{i\vartheta})$ has mean $VR^2$. Fourier coefficient extraction gives $|c_k|R^k\le2VR^2$. The coefficient normalization, the factor two, and the identification $c_k=\kappa_k/k!$ are correct. Hence

$$
|\kappa_k(X_h)|\le2k!V(4W)^{k-2}\qquad(k\ge3).
$$

For $\delta=W/\sqrt V$, the exact summed remainder is

$$
\left|K\!\left(\frac{is}{\sqrt V}\right)+\frac{s^2}{2}\right|
\le\frac{8\delta|s|^3}{1-4\delta|s|}.
$$

At the cutoff $T=(36\delta)^{-1}$, the denominator is at least $8/9$, and the right side is at most $9\delta|s|^3\le s^2/4$. Thus

$$
|\phi(s)-e^{-s^2/2}|\le9\delta|s|^3e^{-s^2/4}
\qquad(|s|\le T).
$$

The exact normal-reference smoothing inequality used in the note was checked in equation (2.16) of [Lebowitz–Pittel–Ruelle–Speer, arXiv version 2](https://arxiv.org/pdf/1408.4153v2). That equation has the integral coefficient $1/\pi$, integration over $[-T,T]$, and error $24/(\pi\sqrt{2\pi}T)$, as quoted. The same source attributes the smoothing inequality to Feller, Section XVI.4, Lemma 2. This audit checked that precise external inequality, not a general normal-approximation theorem from that paper.

Using $\int_{\mathbb R}s^2e^{-s^2/4}\,ds=4\sqrt\pi$ gives

$$
C_*=\frac{36}{\sqrt\pi}
 +\frac{864}{\pi\sqrt{2\pi}}
=130.027827902\ldots<131.
$$

The theorem uses the simpler, conservative constant **150**. This is consistent throughout the reviewed abstract, theorem statements, proof, and scope section. The proof also takes the minimum with one, so no separate small-variance restriction is hidden in the Fourier argument. Neither 150 nor the smaller calculated constant is claimed to be optimal.

## 7. Sharpness is a statement about order

For $m$ disjoint batches of size $B$, each selected with probability $1/m$, time $t=m\log2$ makes their discovery indicators independent Bernoulli variables of parameter $1/2$. Hence

$$
X\stackrel d=B\,\mathrm{Bin}(m,1/2),\qquad V=B^2m/4.
$$

A jump of size $p$ in a distribution function forces Kolmogorov distance at least $p/2$ from every continuous distribution. The maximal binomial atom is asymptotic to $\sqrt{2/(\pi m)}$. Therefore

$$
\liminf_{m\to\infty}\frac{\sqrt V}{B}
d_{\mathrm K}\!\left(\frac{X-\mathbb EX}{\sqrt V},\mathcal N(0,1)\right)
\ge\frac1{2\sqrt{2\pi}}.
$$

This verifies sharpness of the scale $B/\sqrt V$, not of the numerical upper-bound constant. Weighted disjoint singletons of weight $W$ give the same obstruction for the weighted theorem. The example also explains why the result should not be presented as a span-one local limit theorem.

## 8. A rejected Stein-factor shortcut

An alternative Poisson Stein calculation gives the exact positive factor

$$
T=\sum_{i,j}c_{ij}Y_i g_{a_j}(N_j),\quad
g_a(n)=\int_0^1u^n e^{-au}\,du,\quad \mathbb ET=V,
$$

where $a_i$ are species incidence means and $c_{ij}$ are pair-overlap means. The new proof does **not** rely on a bound for $\operatorname{Var}T$.

In particular, the tempting shortcut

$$
\sum_{i,j}c_{ij}(a_i+a_j)
\mathbb E[Y_i g_{a_j}(N_j)^2]\le CV
$$

is false uniformly in the model. To see this, take $M$ high-rate species, each with private singleton intensity $a-1$, and one low-rate species. A common batch containing all $M+1$ species has intensity one. All high incidence means are $a$, the low mean is one, and distinct pair-overlap means are one. Choose $M\sim e^a$. Then

$$
V=M(e^{-a}-e^{-2a})+M(M-1)e^{-2a}(e-1)
  +2Me^{-a}(1-e^{-1})+e^{-1}(1-e^{-1})
$$

stays bounded. The terms with a high first species and low second species alone contribute
$M(a+1)e^{-a}(1-e^{-1})^2\to\infty$ to the proposed left side. This rejects that intermediate estimate only. It is not a counterexample to the optimal-rate theorem, and no claim about the unresolved Stein-factor variance estimate is needed for the analytic proof.

## Scope of this audit

The audit checked the logical and algebraic proof chain, its finite-to-countable passage, its weighted hypotheses, the imported actual-variance argument, the cited smoothing inequality, and the order-sharpness example. Numerical evaluation of the displayed constant supplements these analytic checks; it does not serve as proof of the theorem. Build checks, PDF rendering, finite-model verification scripts, and any subsequent external review are separate validation activities.

## Repository integration checks

The subsequent integration run of `make all` passed all three verification suites and both manuscript validators. The foundation source and PDF remain unchanged. The new companion PDF contains nine pages; all nine were rendered and reviewed for layout, with the complete zero-free induction and Fourier argument also inspected at full-page resolution. No overfull boxes, unresolved references or citations, missing expected page headers, or text-margin violations were found.

The new deterministic suite evaluates nine models, including four weighted models, at three times using 80-decimal-digit arithmetic. It checks 162 cumulant inequalities, 162 independent-set/absence-polynomial identities, 108 weighted Laplace identities, 1,044 induced subgraphs (990 proper), and 8,166 occupation-ratio recursion identities. Its largest reported mean, variance, or partition identity discrepancy is below $1.3\times10^{-80}$; these are high-precision finite calculations, not formal interval certificates. Exact model parameters, tolerances, and script provenance are recorded in [the verification report](../reports/optimal_rate_verification.json).

The [document validation report](../reports/optimal_manuscript_validation.json) records source and PDF fingerprints, label and citation checks, compilation-log checks, and PDF text geometry. The build workflow now compiles and checks both manuscripts. External CI results and the review/merge record are available from the repository's pull request and Actions history.
