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
  itself is rigorous); plus the numerically observed *full-ODE*
  collision-time ratio η ≈ 0.46 (approximately stable, CV ≈ 2.7%, with a weak
  gap-structure dependence), the square-root Speiser scaling, and
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

## Self-consistency: definition and curvature dichotomy

The self-consistency relation does not yield a *universal* upper bound on h₀:
in the regime f″(0) ≥ 0 it gives a **lower** bound, not an upper bound. The
self-consistent depth h_thr = sqrt(2/S_on) is therefore stated as a definition
(the marginal-stability level set f″(0) = 0), with a curvature dichotomy:

- The self-consistency structure consists of
  - a **Definition** of the self-consistent depth
    `h_thr := sqrt(2 / S_on)` (equivalently the level set {f″(0) = 0}), and
  - a **Curvature Dichotomy Proposition**: f″(0) < 0 ⇒ upper bound h₀ < h_thr;
    f″(0) ≥ 0 ⇒ lower bound only.
- The bound **Λ ≤ 0.081** uses the critical-strip constraint h₀ < ½, not the
  self-consistency bound.
- The **structural invariant `h_thr² · S_on = 2`** is exact by
  construction (it *defines* h_thr).

**Layer separation.** Λ ≤ 0.081 (Layer A) rests on the two-zero integral η₂ ≈ 0.487
(rigorous quadrature) and the critical-strip ceiling h₀ < ½; it does not use the
full-ODE ratio η ≈ 0.46 (Layer B). The self-consistent depth `h_thr` is a
**definition** (the marginal-stability level set f″(0) = 0), not a universal upper
bound on off-line zeros.

---

## Scripts

> Filenames below follow the `verify_*.py` convention used in the manuscript;
> adjust to your local names. Each script prints the quantities it verifies.

| Script | Verifies | Layer |
|---|---|---|
| `verify_01_fundamental_identity.py` | V′(h,t) > 0 across the critical strip | A |
| `verify_02_curvature_test.py` | Hadamard residual T₃₀(t), tail accounting | A |
| `verify_03_statistical_balance.py` | Sigma-equivalent of the needed fluctuation | A |
| `verify_04_inner_function.py` | \|Φ_off\| = 1 on the critical line (machine precision) | A |
| `verify_05_concavity.py` | Concavity ∂ₜ²log\|ζ\| < 0, one turning point per gap | A |
| `verify_06_tunneling.py` | Mh₀ = π s/√Cₙ (T-independent); Cₙ ∈ [1.002, 2.306] | B |
| `verify_07_self_consistency.py` | `h_thr² · S_on = 2` to machine precision (exact by definition) | A |
| `verify_08_speiser.py` | Square-root Speiser depth scaling, C̄ ≈ 1.18 | B |
| `verify_09_collision.py` | Full-ODE collision ratio η ≈ 0.46 (CV ≈ 2.7%, weak gap-dependence) | B |
| `verify_10_energy_budget.py` | Λ ≤ 0.081 (two-zero integral) + Poisson budget B(T), E/B → ∞ | A/B |

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
python run_all.py                    # runs all ten verification scripts
python verify_10_energy_budget.py   # Layer A: Lambda <= 0.081; Layer B: E/B -> infinity
python verify_07_self_consistency.py # Layer A: h_thr^2 * S_on = 2
python verify_09_collision.py       # Layer B: full-ODE eta ~ 0.46 (CV ~ 2.7%)
python verify_06_tunneling.py       # Layer B: Mh0 = pi*s/sqrt(Cn), Cn range
```

---

## Citation

Preprint / DOI: see the manuscript. Code DOI on Zenodo (if applicable).

## License

MIT (or as stated in `LICENSE`).
