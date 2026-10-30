# Model, analysis and exact-arithmetic guarantees

This note states the mathematical scope. It does not claim that floating-point implementation is automatically certified, or that the classical spectral inequalities are new.

## 1. Model

There are two equally sized agent populations. A fixed periodic schedule enables one direction at a time, with a positive idle duration after each pulse. Within a phase,

\[
\dot x=B_kx,\qquad B_k=-D+\beta A_k,
\]

where D is positive diagonal, beta is nonnegative and A_k is nonnegative. A_k is either an upper or lower cross-population block; within-population propagation is excluded. The perturbation x is a modelled signal, not a proof that a natural-language decision is correct.

A phase propagator E_k is nonnegative because B_k is Metzler. The period propagator is Phi=E_K ... E_1, with period P. For the fixed periodically repeated schedule, asymptotic stability is equivalent to rho(Phi)<1. Equality is a boundary case. The growth rate is log(rho(Phi))/P. This is not a guarantee for arbitrary schedules or nonperiodic changes.

For diagonal losses d_i and an active cross-edge weight w_ij, the pulse response is

\[
(E_{\rm pulse})_{ij}=\beta w_{ij}\,h\,e^{-d_i h}\,\operatorname{exprel}((d_i-d_j)h),
\]

where exprel(z)=(exp(z)-1)/z and exprel(0)=1. The following idle phase left-multiplies by diag(exp(-d_i h)). Implemented phase operators have diagonal exp(-2d_i h) and one off-diagonal block. An independent implementation using general matrix exponentials validates the formula.

## 2. Baseline bound and warm reuse

For any entrywise positive v and nonnegative Phi, define

\[
\ell(v)=\min_i (\Phi v)_i/v_i,\qquad
u(v)=\max_i (\Phi v)_i/v_i.
\]

Then ell(v) <= rho(Phi) <= u(v). This is the classical Collatz–Wielandt bound. One way to see the upper bound is to use the weighted infinity norm induced by v. The lower bound follows from Phi v >= ell v and iteration (or a positive left Perron vector when irreducible). Positivity of the test vector is necessary; irreducibility is not required for bound validity.

Thus u<1 implies stability, ell>1 implies instability, and otherwise the screen returns unknown. Repeating v <- Phi v / ||Phi v||_infinity refines the bound; intersections of all previously valid intervals remain valid. For primitive Phi, normalized power iteration converges to a Perron vector. No finite uniform pass bound is claimed for arbitrary positive systems.

Cold starts from ones; warm starts from a positive vector retained from the previous version. The retained vector may be approximate: its approximation error affects tightness, not exact-arithmetic validity of the inequality. Floating-point errors in evaluating the action and bounds are a separate issue.

With a structured block phase, one matrix-free period action costs O(K n²) for dense cross-kernels; a budget of b passes costs O(b K n²). Dense product construction costs O(K n³), plus eigendecomposition. Matrix-free cold is therefore an essential baseline: it separates the benefit of structure from the benefit of previous-version reuse.

## 3. Exploratory phase-cache delta screen

Let T_k be the old phase operators and T'_k the candidate phase operators. Cache the old trajectory of any v_0>0:

\[
v_k=T_kv_{k-1},\quad k=1,\dots,K.
\]

The cache contains K+1 positive vectors. Compute the old residual ratios

\[
q_- = \min_i(v_K)_i/(v_0)_i,\qquad
q_+ = \max_i(v_K)_i/(v_0)_i.
\]

For each changed phase, compute a local comparison using only its delta:

\[
r_{ki}=1+\frac{[(T'_k-T_k)v_{k-1}]_i}{(v_k)_i},\qquad
\alpha_k=\min_i r_{ki},\quad \omega_k=\max_i r_{ki}.
\]

Because r represents T'_k v_{k-1}/v_k, it is positive under the assumed positive diagonal phase model. We have

\[
\alpha_k v_k\le T'_k v_{k-1}\le\omega_kv_k.
\]

**Composition bound (exact arithmetic).** Positivity and induction give

\[
\left(\prod_k\alpha_k\right)v_K
\le\Phi'v_0\le
\left(\prod_k\omega_k\right)v_K.
\]

Consequently,

\[
q_-\prod_k\alpha_k\le\rho(\Phi')\le q_+\prod_k\omega_k.
\]

**Proof.** Suppose after k-1 candidate phases the propagated vector lies between the products of previous local lower/upper factors times v_{k-1}. Apply the nonnegative T'_k and then the local comparison inequality. This proves the vector inequality inductively. Divide by v_0 componentwise and apply the spectral bound from Section 2. No equality v_K=q v_0 is required.

For an update affecting r rows of a cross-kernel and an unchanged schedule, delta-vector products use O(K r n) arithmetic. The local extrema still scan O(K n) entries. The cache costs O(K n) extra storage, in addition to existing kernels. In the current code, preparing candidate kernels and extracting delta rows still cost O(n²). Therefore the end-to-end implementation is not O(Krn) without qualifications. Direction changes can require full block-vector products and have weaker locality.

The local min/max factors may combine worst cases that are not jointly attainable. The cache bound is therefore potentially much more conservative than evaluating Phi' v_0 directly. The hybrid first attempts this screen and falls back to warm matrix-free CW when unknown. Failure to conclude is not acceptance.

This construction uses established positive-operator comparison. Its novelty relative to existing compositional and incremental stability methods must be checked before a publication claim. The experiment introduces it as an exploratory candidate, after the main CW study; it is not a preregistered confirmatory comparison.

## 4. Numerical and validation scope

All main and exploratory screens use float64 and a 1e-10 guard around rho=1. Such a guard is not a general rounding-error proof. Dense reference disagreement checks are numerical evidence, not formal verification.

Thirty small-system candidates are independently re-evaluated using 60-digit interval arithmetic and analytic phase response. Each candidate uses a strictly positive vector after sixteen floating-point refinement passes, then the full phase action is independently recomputed as an interval. This produces interval-based decisions for those supplied binary64 input values considered exact. It does not cover data uncertainty or automatically validate the cache layer.

The cache implementation is cross-checked against a scalar implementation on all 3,600 main candidates and against the dense spectral reference. Scalar timing outputs are retained separately after the vectorized implementation was introduced, to expose the implementation revision.

## 5. Publication positioning

The present artifact supports a paper about version-aware feedback-stability analysis within a restricted periodic positive-system model. It cannot support claims about arbitrary MAS evolution, all schedules, correctness of LLM outputs, automated repair, or safety of production systems. A paper needs either a distinct, literature-validated analysis contribution, or an executable application case study with a justified model-to-system mapping, in addition to these synthetic results.

Foundational background: Perron–Frobenius theory, Collatz–Wielandt inequalities, positive switched systems, and periodic linear-system stability. The supplied supplement motivates chronology-sensitive analysis; this artifact does not reproduce its compressed-log certificate algorithm.
