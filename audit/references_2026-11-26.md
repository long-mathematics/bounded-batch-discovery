# Reference and Janisch application audit

Audit date: 26 November 2026. This is a source and theorem-statement audit, not
a formal verification of either manuscript and not an exhaustive priority search.
No repository source files were changed by this audit.

## Conclusions

The inspected primary sources support the reference claims in the current
manuscript and the bounded-degree corollary in `extensions/berry_esseen.md`.
No substantive mathematical correction to that corollary was identified.
The main maintenance changes recommended here are permanent source links,
GitHub-compatible mathematical delimiters, and a clear separation between
historical and newly rerun validation records.

In particular, **do not change the Haeusler parameter range from p >= 1 to
p > 1**: Mourrat's Theorem 1.1 explicitly includes p = 1, and the manuscript
correctly says its applications use fixed p >= 2.

## 1. Janisch source and exact constants

The public repository was cloned and inspected at commit
`8a2c9e9a4eee068413f0a5fdfbf8e16172870c8a`. Its commit timestamp is
26 November 2026, 19:23:22 +01:00. The LaTeX manuscript date is
26 November 2026.

Permanent inspected sources:

- [Manuscript source](https://github.com/maximilian-janisch/berry-esseen-dependency-graphs/blob/8a2c9e9a4eee068413f0a5fdfbf8e16172870c8a/paper/janisch-lehericy-proof.tex).
- [Compiled manuscript](https://github.com/maximilian-janisch/berry-esseen-dependency-graphs/blob/8a2c9e9a4eee068413f0a5fdfbf8e16172870c8a/paper/janisch-lehericy-proof.pdf).
- [Formal assumptions and conclusions](https://github.com/maximilian-janisch/berry-esseen-dependency-graphs/blob/8a2c9e9a4eee068413f0a5fdfbf8e16172870c8a/JanischLehericy/Statement.lean).
- [Normalized theorem](https://github.com/maximilian-janisch/berry-esseen-dependency-graphs/blob/8a2c9e9a4eee068413f0a5fdfbf8e16172870c8a/JanischLehericy/Theorem.lean).
- [Published-form theorem](https://github.com/maximilian-janisch/berry-esseen-dependency-graphs/blob/8a2c9e9a4eee068413f0a5fdfbf8e16172870c8a/JanischLehericy/Published.lean).
- [Upstream verification instructions](https://github.com/maximilian-janisch/berry-esseen-dependency-graphs/blob/8a2c9e9a4eee068413f0a5fdfbf8e16172870c8a/VERIFY.md).

Theorem 1 assumes a finite family of centered, square-integrable variables,
variance one for their sum, and a genuine dependency graph: disjoint vertex
sets with no joining edge have independent random vectors. For a positive
integer q bounding all closed neighborhood sizes, it defines

$$
\ell(u)=u^2\min(1,u),\qquad
\beta=q^{-1}\sum_i\mathbb E\ell(q|X_i|)
$$

and states the Kolmogorov bound $\min(1,250000\beta)$.
Corollary 2 gives, for $2<\delta\le3$,

$$
\min\left(1,250000q^{\delta-1}\sum_i\mathbb E|X_i|^\delta\right).
$$

Theorem 3 uses uncentered/raw moments and constant 2,000,000. The extra
factor eight is the centering estimate. Our discovery application already
uses centered Bernoulli variables, so its use of 250,000 is correct.

The inspected Lean definitions use the CDF of `gaussianReal 0 1`, ordinary
variance, measurable square-integrable centered summands, and independent
separated blocks. The hypotheses and conclusions inspected at the theorem
boundary match the mathematical statement relevant here. The repository pins
Lean 4.24.0 and mathlib v4.24.0. No Lean build, axiom audit, or kernel replay
was run as part of this audit. Upstream reports about those checks must remain
attributed to upstream. A source inspection is not a substitute for replay.

## 2. Verification of the bounded-degree discovery corollary

The reasoning in note lines 22–38 is sound, conditional on Janisch's stated
theorem. For each supported finite batch C, Poisson splitting supplies an
independent batch-type process. Each discovery indicator uses exactly the
processes whose batches contain its species. If no edge joins two sets of
species, no positive-rate batch meets both sets. Their indicator vectors are
therefore functions of disjoint independent families. This verifies graph
dependence, rather than merely pairwise independence.

For a finite species restriction J, put
$X_i=(I_i-p_i)/\sqrt{V_J}$ and $q=D+1$. Applying the centered third-moment
corollary gives the displayed expression. The elementary identity

$$
\mathbb E|I_i-p_i|^3
=p_i(1-p_i)\bigl(p_i^2+(1-p_i)^2\bigr)
\le p_i(1-p_i)
$$

is correct. The joint-absence formula makes all covariances nonnegative,
so $\sum_{i\in J}p_i(1-p_i)\le V_J$. The resulting bound is
$\min\{1,250000(D+1)^2/\sqrt{V_J}\}$.

For an increasing exhaustion of countably many species, counts converge
almost surely and in $L^2$: every restricted count and the full count are
bounded by $BN(t)$, which has a finite second moment. Consequently means
and variances converge. Eventually the restricted variances are positive
whenever the limiting V is positive. The standardized sums converge in law,
and the uniform inequalities against the continuous normal CDF pass to the
limit, including possible atoms, by taking nearby thresholds. Thus the
countable extension does not need a uniform bound on the number of species.

The degree/batch-size distinction is also correct. Uniform batches of two
distinct species on n species have degree n-1, marginal rate 2/n, and
pair-to-marginal ratio 1/(n-1). The identity
$\sum_{j\ne i}\lambda_{ij}\le(B-1)\lambda_i$ controls total overlap
weight; it does not bound degree.

The note correctly leaves the universal $CB/\sqrt V$ estimate open here.
Its lattice-order argument can be made completely concrete: choose N
disjoint packages of B species uniformly and use time t=cN. Then
$X=B\operatorname{Bin}(N,1-e^{-c})$, with
$V=B^2Np(1-p)$. The largest atom has order $N^{-1/2}$ and any continuous
comparison has Kolmogorov error at least half that atom. This yields the
claimed necessary order $B/\sqrt V$. This does not prove the proposed upper
bound or a novelty claim.

### Recommended note edits

1. Replace the moving Janisch repository-only citation with the pinned
   manuscript/commit above, identifying Theorem 1 and Corollary 2.
2. Keep the explicit statement that the Lean formalization has not been
   independently replayed for this project.
3. Convert `\(...\)` and `\[...\]` to GitHub's `$...$` and `$$...$$` in
   Markdown prose. This affects presentation, not the mathematics.
4. If the current paper is revised during repository setup, distinguish the
   preserved November 12 version from the current revision in the opening
   paragraph and closing provenance sentence.

## 3. Haeusler/Mourrat

Inspected [Mourrat, arXiv:1103.5050v2](https://arxiv.org/html/1103.5050v2),
Theorem 1.1, and downloaded the corresponding PDF. The theorem states
p in [1,infinity), no bounded-increment assumption, and

$$
D(\mathbf X)\le C_p\left(
\|V^2(\mathbf X)-1\|_p^p
+s^{-2p}(\mathbf X)\sum_i\|X_i\|_{2p}^{2p}
\right)^{1/(2p+1)}.
$$

This is precisely the external inequality used at current manuscript lines
588–601. Haeusler is credited through the explicitly inspected Mourrat
formulation. We did not independently check Haeusler's original proof.

## 4. Janson and Eriksson

### Janson

Inspected the [published article](https://pmc.ncbi.nlm.nih.gov/articles/PMC6405020/),
Problem 6.9. It asks for a vertex-count CLT under divergence of that count's
variance. The text immediately following starts the separate discussion of
dust and attached stars. The manuscript's restriction to a fixed law on a
countable species space is appropriate. This source audit verifies what the
problem asks, not the correctness of our proposed solution.

### Eriksson graph paper

Inspected [arXiv:2501.09511v2](https://arxiv.org/html/2501.09511v2) and downloaded
its PDF. The header verifies 6 May 2026. Theorem 14 is the vertex CLT under
eventual forever connectedness; Proposition 16 handles a support consisting
of isolated edges; Proposition 17 compares the vertex and independent-urn
variances. Appendix A, Lemma 30 proves $f(t)^2/f(2t)\to0$ under the stated
infinite-support setting and credits the earlier urn argument. The current
manuscript's theorem numbers and attributions agree with this pinned version.

### Eriksson unseen-species paper

Downloaded and inspected [arXiv:2602.08769v4](https://arxiv.org/pdf/2602.08769v4).
The PDF header and [arXiv version record](https://arxiv.org/abs/2602.08769)
verify 9 September 2026. Lemma 29, printed page 22, gives concentration for
the count of species seen at least r times, assuming bounded batch arity B.
Its scale contains $((B-1)r+1)$ and the expected threshold count. The proof
on printed pages 71–72 constructs a bounded size-bias coupling and invokes
the stated prior coupling theorem. This supports the manuscript's warning
against claiming the first incidence-data concentration result. We did not
independently re-prove that external coupling theorem.

## 5. Relation to the recovered November 12 reference audit

Read the recovered `reference_audit.md` from the revised November 12 bundle.
It already distinguishes exact theorem checks, bibliographic checks,
contextual citations, and access limitations. Its Haeusler, Janson, and
Eriksson entries are corroborated by the current inspection above.

Other original entries discuss Boucheron–Lugosi–Massart, Finner/GLSS,
Esary–Proschan–Walkup, Dembo–Zeitouni, Dembo, Dutko, and Good–Toulmin.
Their full-text checks are historical audit statements; this follow-up did
not rerun every one. Preserve that distinction in the public audit record.
In particular, do not turn the original explicit limitations for Finner's
full text or the subscription-gated book theorem into a claim of direct
inspection today.

No external PDFs or Lean sources need to be copied into our repository.
The pinned links above suffice for reproducibility of the statement audit.
