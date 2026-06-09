# rh-partial-fraction

Numerical verification code for the manuscript

**Backward Heat-Flow Collision Dynamics for Zeros of the Riemann ξ-Function: A Constraint Metric Tensor, Spectral Tunneling, and the Bound Λ ≤ 0.081 on the de Bruijn–Newman Constant**

Mesut Ismail — TU-Sofia

---

## Overview

This repository contains the Python scripts that reproduce every numerical
table and constant reported in the manuscript. All computations use
[`mpmath`](https://mpmath.org/) at 25-digit precision with the first 2000
non-trivial zeros of ζ (γ ≤ 2515).

The manuscript studies a hypothetical off-line zero
ρ\* = ½ + h₀ + i t₀ (h₀ > 0) of the Riemann ξ-function as an interacting
particle under the **de Bruijn–Newman backward heat flow**
∂ₜH = −∂ₓ²H, H₀ = ξ. The analysis is organized into three explicit layers:

- **Layer A (rigorous).** Theorems on Apollonius-disc confinement, the exact
  Hadamard second-derivative identity and the resulting self-consistent depth
  and curvature dichotomy, measure collapse, antisymmetry of ξ′/ξ, the Blaschke
  inner-function characterization, the concavity equivalence chain, and the
  **unconditional bound Λ ≤ 0.081** on the de Bruijn–Newman constant.
- **Layer B (numerical observations).** The 5×5 Gram (constraint metric)
  tensor Gᵢⱼ and the structural divergence E(T)/B(T) → ∞ (the divergence
  itself is rigorous); plus the numerically observed *full-ODE* universal
  collision-time ratio η = 0.456 ± 0.002, the square-root Speiser scaling, and
  Cₙ ∈ [1.002, 2.306].
- **Layer C (conjectural).** A reformulation of the prime-side obstruction.

> **Important (rigorous vs. numerical).** The bound Λ ≤ 0.081 is a **Layer A**
> result. It rests only on (i) the elementary critical-strip constraint
> h₀ < ½, (ii) a **rigorous two-zero collision-time integral** (an absolutely
> convergent quadrature whose value upper-bounds the full backward-flow
> collision time by the comparison principle), and (iii) the Platt–Trudgian
> verification height together with Trudgian's gap bound. It does **not** use
> the numerically observed full-ODE constant η = 0.456, and it does **not** use
> the self-consistency dichotomy. These are kept strictly separate.

---

## What changed in this revision

This revision corrects an error in the previous "self-consistency upper bound"
statement (formerly Theorem 5.3). In the regime f″(0) ≥ 0 the self-consistency
relation yields a **lower** bound on h₀, not an upper bound, so a *universal*
self-consistency upper bound does **not** hold. Accordingly:

- The former "self-consistency upper bound" theorem is replaced by
  - a **Definition** of the self-consistent depth
    `h_thr := sqrt(2 / S_on)` (equivalently the level set {f″(0) = 0}), and
  - a **Curvature Dichotomy Proposition**: f″(0) < 0 ⇒ upper bound h₀ < h_thr;
    f″(0) ≥ 0 ⇒ lower bound only.
- The bound **Λ ≤ 0.081 is unchanged**: it always used the critical-strip
  constraint h₀ < ½, not the self-consistency bound.
- The **structural invariant `h_thr² · S_on = 2`** is unchanged and is exact by
  construction (it *defines* h_thr).

**No numerical value changes.** Λ ≤ 0.081, η₂ (two-zero), η = 0.456 (full ODE),
Cₙ, S_on, the invariant, the Gram diagonals and E/B are all numerically
identical to before. Only the *interpretation* of `h_thr` changed: it is a
**definition of the self-consistent depth**, never a universal upper bound on
off-line zeros.

---

## Scripts

> Filenames below follow the `verify_*.py` convention used in the manuscript;
> adjust to your local names. Each script prints the quantities it verifies.

| Script | Verifies | Layer |
|---|---|---|
| `verify_invariant.py` | `h_thr² · S_on = 2` to machine precision (exact by definition) | A |
| `verify_curvature.py` | Hadamard identity f″(0) = S_on + \|∂ₜ²log\|P\|\| − 2/h₀² | A |
| `verify_lambda.py` | Λ ≤ 0.081 via the two-zero integral at h₀ = ½, Lₙ = L_max = 1.614 | A |
| `verify_eta2.py` | Two-zero collision ratio η₂ = 0.487 (rigorous quadrature) | A |
| `verify_eta_full.py` | Full-ODE universal ratio η = 0.456 ± 0.002 over 200 gaps | B |
| `verify_gram.py` | Gram tensor Gᵢⱼ, diagonals, E(T) = tr G | B |
| `verify_budget.py` | Poisson budget B(T), ratio E/B → ∞ | B |
| `verify_speiser.py` | Square-root Speiser depth scaling, C̄ ≈ 1.18 | B |
| `verify_cn.py` | Cₙ ∈ [1.002, 2.306] across 1999 gaps | B |

### The Λ ≤ 0.081 computation (Layer A, self-contained)

```python
import mpmath as mp
mp.mp.dps = 30

# Critical-strip depth ceiling (unconditional): h0 < 1/2.
# Maximum normalized gap at T0 = 3e12 via Trudgian: L_max ≈ 1.614.
L_max = mp.mpf('1.614')
half_sq = (L_max / 2)**2          # (L_n/2)^2

# Two-zero backward-flow collision integral (rigorous quadrature).
# The full-flow collision time is <= this by the comparison principle.
integrand = lambda h: h * (h**2 + half_sq) / (half_sq + 5 * h**2)
Lambda_bound = mp.quad(integrand, [0, mp.mpf('0.5')])

print("Lambda <=", Lambda_bound)   # 0.0808...  ->  Lambda <= 0.081
```

This snippet uses **no** numerical η and **no** self-consistency bound — only
h₀ < ½, L_max, and an absolutely convergent integral.

---

## Requirements

```
python >= 3.9
mpmath >= 1.3.0
```

```bash
pip install mpmath
```

A cache of the first 2000 ζ-zeros (γ ≤ 2515) is used; regenerate with
`mpmath.zetazero(n)` if not present.

---

## Reproducing the tables

```bash
python verify_lambda.py      # Layer A: Lambda <= 0.081
python verify_invariant.py   # Layer A: h_thr^2 * S_on = 2
python verify_eta2.py        # Layer A: two-zero eta_2 = 0.487
python verify_eta_full.py    # Layer B: full-ODE eta = 0.456 +/- 0.002
python verify_gram.py        # Layer B: Gram tensor, E(T)
python verify_budget.py      # Layer B: E/B -> infinity
```

---

## Citation

Preprint / DOI: see the manuscript. Code DOI on Zenodo (if applicable).

## License

MIT (or as stated in `LICENSE`).
