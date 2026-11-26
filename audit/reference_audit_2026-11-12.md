# Reference and dependency audit

**Revision date:** 12 November 2026.

This record distinguishes direct checks of theorem statements from bibliographic
checks and contextual references. The manuscript gives complete proofs of its
new inequalities and reductions; it does not attribute those derivations to
papers that contain only related machinery. This is not an exhaustive priority
search. External PDFs are not redistributed in the bundle.

## External mathematical inputs

### Haeusler (1988), recorded by Mourrat (2013)

- Haeusler: *On the rate of convergence in the central limit theorem for
  martingales with discrete and continuous time*, Annals of Probability 16(1),
  275–299. DOI: `10.1214/aop/1176991901`.
- Mourrat: *On the rate of convergence in the martingale central limit theorem*,
  Bernoulli 19(2), 633–645. DOI: `10.3150/12-BEJ417`.
- Checked the exact statement of **Mourrat, Theorem 1.1**, in the author/journal
  reprint at https://arxiv.org/pdf/1103.5050 (arXiv v2, 21 March 2013).
- The inequality is available for p in [1,infinity), and does **not** impose a
  bounded-increment hypothesis. The manuscript uses p>=2 for its moment bounds.
  It attributes the result to Haeusler and explicitly identifies the inspected
  formulation in Mourrat; it does not claim to have re-audited Haeusler's original
  proof.
- The proof of **Mourrat, Theorem 1.5** was checked for the stopping-and-augmentation
  precedent. It is not invoked as a black box on deterministic-grid increments.
  The manuscript supplies its own Gaussian completion and smoothing proof.

### Finner (1992) and Gavinsky–Lovett–Saks–Srinivasan (2015)

- Finner: *A generalization of Holder's inequality and some probability
  inequalities*, Annals of Probability 20(4), 1893–1901.
  DOI: `10.1214/aop/1176989534`.
- GLSS: *A tail bound for read-k families of functions*, Random Structures &
  Algorithms 47(1), 99–108. DOI: `10.1002/rsa.20532`.
- Checked the read-k formulation and its dependence structure in
  https://arxiv.org/pdf/1205.1478, and publication metadata. The original Finner
  full text was not obtained through the publisher gateway. Its bibliographic
  information and role are corroborated by the read-k literature.
- No uninspected theorem hypothesis is needed: Lemmas 2.1–2.2 prove the particular
  product-space inequalities used, including the countable-coordinate passage.
  These abstract inequalities are presented as classical machinery, not new
  discoveries.

### Boucheron–Lugosi–Massart (2000 and 2009)

- 2000: *A sharp concentration inequality with applications*, Random Structures
  & Algorithms 16(3), 277–292.
  DOI: `10.1002/(SICI)1098-2418(200005)16:3<277::AID-RSA4>3.0.CO;2-1`.
- Checked **Theorem 1** in the author preprint
  https://econ-papers.upf.edu/papers/376.pdf, dated 15 April 1999; the bibliography
  correctly records the 2000 journal publication, not the preprint year.
- 2009: *On concentration of self-bounding functions*, Electronic Journal of
  Probability 14, 1884–1899. DOI: `10.1214/EJP.v14-690`.
- Checked the all-real sub-Poisson moment bound recalled on **p. 1885** and the
  final lower-tail inequality of **Theorem 1**, at (a,b)=(1,0), in the journal PDF:
  https://www.maths.tcd.ie/EMIS/journals/EJP-ECP/article/download/690/690-2335-1-PB.pdf.
- The application is to D_n/B. The endpoint and out-of-range lower-tail events
  are treated explicitly. This is a mean-proxy fixed-sample result, not an
  actual-variance concentration theorem for fixed sampling.

### Esary–Proschan–Walkup (1967)

- *Association of random variables, with applications*, Annals of Mathematical
  Statistics 38(5), 1466–1474. DOI: `10.1214/aoms/1177698701`.
- Checked **Theorem 2.1** (independent variables are associated) and **property P4**
  (nondecreasing functions preserve association), p. 1467, in the original
  article scan at https://www.stat.cmu.edu/~brian/720/week02/esary-proschan-walkup-1967.pdf.
- Used for the covariance lower bound for threshold counts. The manuscript also
  recalls the one-coordinate covariance identity and product induction.

### Dembo–Zeitouni, second edition, corrected printing (2010)

- *Large Deviations Techniques and Applications*, Stochastic Modelling and
  Applied Probability 38. DOI: `10.1007/978-3-642-03311-7`.
- The publisher record confirms that this is the corrected printing of the
  second edition originally printed in 1998:
  https://link.springer.com/book/10.1007/978-3-642-03311-7.
- The cited **Theorem 2.3.6 (Gartner–Ellis)** is the standard finite-dimensional
  cumulant-generating-function theorem; the same theorem locator is explicitly
  referenced in Dembo's 1996 paper, p. 14. The book's subscription-gated theorem
  page was not directly retrieved in this audit.
- The application is the uncomplicated full-domain quadratic case: the manuscript
  proves the scaled log moment-generating function converges to lambda^2/2 for
  every real lambda, including the exceptional-event control under the tilt.
  The limiting function is finite, differentiable and essentially smooth on R,
  with quadratic convex conjugate. No boundary-of-domain hypothesis is being
  omitted.

### Dembo (1996)

- *Moderate deviations for martingales with bounded jumps*, Electronic
  Communications in Probability 1, 11–17. DOI: `10.1214/ECP.v1-973`.
- Checked **Proposition 1 and Corollary 1(a)**, p. 13, in the journal PDF:
  https://www.maths.tcd.ie/EMIS/journals/EJP-ECP/article/download/973/973-3311-1-PB.pdf.
- Cited as the established bounded-jump/exponential-bracket framework. Its
  process-level hypotheses are not silently imported into the manuscript's
  triangular-array theorem: that theorem has its own exponential-tilting proof.

## Model history and comparison

### Janson (2018)

- *On edge exchangeable random graphs*, Journal of Statistical Physics 173,
  448–484. DOI: `10.1007/s10955-017-1832-9`.
- Checked **Problem 6.9** in https://arxiv.org/pdf/1702.06396. It asks for
  asymptotic normality of the discovered-vertex count. The paper's separate
  discussion of more general dust mechanisms is not included in our theorem's
  scope. The result in this manuscript is for a fixed law on a countable species
  space.

### Eriksson, graph paper, version 2

- *Edge exchangeable graphs: connectedness, Gaussianity and completeness*,
  arXiv:2501.09511v2, **6 May 2026**.
- Checked the user-supplied paper against https://arxiv.org/pdf/2501.09511v2.
- **Theorem 14:** vertex CLT in the eventually forever connected regime.
- **Proposition 16:** isolated-edge case.
- **Proposition 17 / Corollary 18:** covariance comparison and connected-regime
  variance equivalence.
- **Lemma 30 (Appendix A):** f(t)^2/f(2t) tends to zero.
- The draft credits this existing coupling and de-Poissonization work, rather
  than presenting its exponential clocks or derivative estimate as new.

### Dutko (1989)

- *Central limit theorems for infinite urn models*, Annals of Probability 17(3),
  1255–1263. DOI: `10.1214/aop/1176991268`.
- Bibliographic record checked. The derivative estimate is credited as the
  classical urn argument and is directly checked through Eriksson's Lemma 30,
  where the earlier attribution is explicit. The manuscript repeats the proof;
  it does not assign a guessed theorem number in Dutko's original article.

### Eriksson, unseen-species paper, version 4

- *The unseen species problem revisited*, arXiv:2602.08769v4,
  **9 September 2026**. Source: https://arxiv.org/pdf/2602.08769v4.
- Checked **Lemma 29**, concerning concentration for multiple-coverage counts.
  The comparison in the draft acknowledges that mean-based concentration already
  exists. It does not claim the first concentration result for incidence data,
  or infer unknown-law prediction calibration from an oracle CLT.
- These are citations to specified versions, not an assertion that no newer
  paper or revision exists.

### Good–Toulmin (1956)

- *The number of new species, and the increase in population coverage, when a
  sample is increased*, Biometrika 43(1–2), 45–63.
  DOI: `10.1093/biomet/43.1-2.45`.
- Publication title, authors and metadata checked against the publisher record:
  https://academic.oup.com/biomet/article-abstract/43/1-2/45/334874.
- Included for the named class of truncated occupancy estimators in the scope
  discussion. No prediction-error theorem or optimality claim is attributed to
  this citation in the present manuscript.

## Conclusions of the citation audit

The bibliography contains 14 entries, all cited. Exact external theorem
formulations used in the martingale and self-bounding steps have been checked.
The proof of the Gaussian-completion refinement is supplied in the draft, not
attributed to a nonexistent theorem in the literature. Access limitations for
some original historical sources are recorded above instead of being concealed.
