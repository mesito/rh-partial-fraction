# RH Partial-Fraction Numerical Verifications

Code for reproducing numerical results in:

**"Arithmetic Non-Realizability of Off-Line Blaschke Factors
in the Riemann Zeta Function: Constraint Metric Tensor and the Bound Λ ≤ 0.081"** — Version 13, Mesut Ismail, May 2026

Submitted to *Journal of Number Theory* (MS: JNTH-D-26-00613, 12 May 2026).

## Requirements

- Python 3.8+
- mpmath ≥ 1.3: `pip install mpmath`

## Usage

```bash
python run_all.py             # Full verification (~45 min with 2000 zeros)
python run_all.py --quick     # Quick mode (~10 min, 200 zeros)

# Individual verifications:
python verify_01_fundamental_identity.py
python verify_02_curvature_test.py
# ... etc.
```

## Architecture

Zeros are loaded in batches of 200 via `mpmath.zetazero(n)` at 25-digit
precision and cached to `zeros_cache.pkl`. The default loads 2000 zeros
(10 batches of 200); adjust `N_BATCHES` in `config.py` for faster/slower runs.

All computations use `mpmath.mp.dps = 25` (25 significant decimal digits).

## File Map

| File | Section | What it verifies |
|------|---------|------------------|
| `config.py` | — | Shared config, zero loading, S_on computation |
| `verify_01_*` | §17.1 | V'(h,t) > 0 at 70 grid points (Theorem 4.1) |
| `verify_02_*` | §17.2 | Curvature residual T₃₀, S/N > 10⁴ (Theorem 10.1) |
| `verify_03_*` | §17.3 | Sigma-equivalents ≫ 1 (Proposition 11.1) |
| `verify_04_*` | §17.4 | \|Φ_off\| = 1 on CL exactly (Theorem 12.1) |
| `verify_05_*` | §17.5 | Concavity d²/dt² < 0, 1 turning point/gap (Theorem 14.1) |
| `verify_06_*` | §17.6 | Tunneling Mh₀ = πs/√Cₙ, Cₙ statistics (Theorem 15.1) |
| `verify_07_*` | §17.7 | Self-consistency h₀²·S_on = 2 to machine precision |
| `verify_08_*` | §17.8 | Speiser depth scaling, √δ law (Theorem 13.1) |
| `verify_09_*` | §17.9 | Collision-time ODE, universal η ≈ 0.456 (Theorem 16.1) |
| `verify_10_*` | §§18–27 | Gram tensor E(T), budget B(T), E/B → ∞ (Theorem 21.2) |
| `verify_paper2.py` | — | 50-digit verification for Paper #2 (Λ ≤ 0.078) |
| `run_all.py` | — | Master runner |

## Verified Results (2000 zeros, γ ≤ 2515)

| Verification | Result (2000 zeros) | Paper claim |
|---|---|---|
| V' > 0 | All 70 points ✓ | All 70 points |
| Curvature S/N | > 13,000 ✓ | > 10⁴ |
| \|Φ_off\| on CL | 1.000000000000000 ✓ | 1 exact |
| Concavity | All 50 gaps concave, 1 TP each ✓ | All gaps |
| Cₙ range | [1.002, 2.306], mean 1.29 ✓ | [1.002, 2.306], mean 1.29 |
| h₀² · S_on | 2.0000000000 ± 0 ✓ | 2.000000 ± 10⁻¹⁴ |
| Mh₀/s | 2.80 ± 0.19 ✓ | 2.80 ± 0.20 |
| Speiser C̄ | 1.04 ± 0.03 | 1.04 ± 0.03 |
| η (collision, full ODE) | 0.456 ± 0.002 | 0.456 ± 0.002 |
| η₂ (two-zero ODE) | 0.487 | 0.487 |
| E/B | grows as log⁵ T ✓ | Theorem 21.2 |
| **Λ (unconditional)** | **≤ 0.081** ✓ | **≤ 0.081** |

## Key Results

### Λ ≤ 0.081 (unconditional)

The main quantitative result (Theorem 22.1). Combines three ingredients:

1. **Platt–Trudgian verification**: all zeros with γ ≤ 3 × 10¹² lie on the critical line.
2. **Self-consistency upper bound** (Theorem 5.3): h₀ ≤ Lₙ/(2√Cₙ) with Cₙ ≥ 1.
3. **Two-zero collision ODE** (Lemma 22.1): η₂ = 0.487, with corrected formula
   dh/dt = −1/h − 4h/(h² + (L/2)²).

The two-zero ODE is a rigorous upper bound on collision time by comparison principle:
the full ODE (with all on-line zeros) has greater acceleration, hence faster collision.
At the critical-strip binding depth h₀ = 1/2 with L_max ≈ 1.614:

```
Λ ≤ τ₂(1/2, 1.614) = ∫₀^{1/2} h(h²+0.651)/(5h²+0.651) dh ≈ 0.081
```

This improves the Platt–Trudgian bound (Λ ≤ 0.20) by a factor of 2.5 via an
independent method — the first new approach to the de Bruijn–Newman constant
since Polymath15 (2019).

### Constraint Network (a)–(k)

Version 13 establishes eleven constraints that any hypothetical off-line zero must
simultaneously satisfy:

| Constraint | Description | Status |
|---|---|---|
| (a) Confinement | Apollonius disc of radius R → 0 | Rigorous |
| (b) Sign barrier | h₀ ≥ c/log T | Conditional |
| (c) Self-consistency | h₀ ≤ Lₙ/(2√Cₙ) ≤ C'/log T | Rigorous |
| (d) Measure collapse | Total disc area → 0 as T^{α−2} | Rigorous |
| (e) Curvature detection | S/N > 10⁴ in Hadamard expansion | Rigorous |
| (f) Statistical rarity | meas{V'<0}/T ≤ exp(−cT²log²T) | Heuristic |
| (g) Invisibility | \|Φ_off\| = 1 on CL | Rigorous |
| (h) Speiser pair | Forced ξ' zero at √δ depth | Rigorous |
| (i) Concavity | d²/dt² log\|ζ\| < 0 in every gap | Rigorous |
| (j) Tunneling parameter | Mh₀ = πs/√Cₙ, T-independent | Rigorous |
| (k) Heat-flow collision | η = 0.456 gives Λ_local ≤ C/log² T | Rigorous/Numerical |

### Structural Inequality E/B → ∞ (Theorem 21.2)

The Speiser-free bound E/B ≥ π log³ T / 2 holds at all depths including
self-consistency, and is independent of the Speiser scaling constant.

## Three-Layer Structure

The paper is organized into three parts with explicit logical status:

- **Layer A (rigorous):** Theorems proved from standard analytic number theory
  (Hadamard factorization, Levinson bound, Titchmarsh bounds, Ingham density estimate).
  Includes Λ ≤ 0.081 unconditional.

- **Layer B (numerical observations):** Gram tensor, collision-time universality
  (η = 0.456 ± 0.002), Speiser scaling (C̄ ≈ 1.04). Structural inequality E/B → ∞
  is rigorous; numerical constants are observations at 25-digit precision.

- **Layer C (conjectural):** Arithmetic Non-Realizability principle. Nontrivial
  off-line Blaschke factors of ξ lie outside the Euler-realizability set R_ζ.
  Supported by Constraint Necessity Lemma and Euler Realizability Bound.

## Changes from v12 → v13

| Item | v12 | v13 | Impact |
|------|-----|-----|--------|
| **Λ bound** | Λ ≤ 0.200 | **Λ ≤ 0.081** | Factor 2.5 improvement |
| **η₂ (two-zero)** | 0.615 | **0.487** | Corrected ODE formula (4h not 4h²) |
| **ODE formula** | dh/dt = −1/h − 4h²/(h²+(L/2)²) | **dh/dt = −1/h − 4h/(h²+(L/2)²)** | Formula correction; numerics matched 0.081 |
| η (full ODE) | 0.459 ± 0.003 | **0.456 ± 0.002** | Refined computation |
| Speiser constant C̄ | 1.18 ± 0.05 | **1.04 ± 0.03** | Refined computation |
| New: Caster effect | — | §8, Theorem 8.1 | Mean drift identity |
| New: Antisymmetry | — | §9, Theorem 9.1 | ξ'/ξ antisymmetry |
| Part III | Arithmetic Closure (v11 style) | **Arithmetic Non-Realizability** | Fully reworked: Euler realizability framework |
| Submission | Preprint only | **Submitted to JNT** | MS: JNTH-D-26-00613 |

## Citation

```bibtex
@article{ismail2026partialfraction,
  title   = {Arithmetic Non-Realizability of Off-Line Blaschke Factors
             in the Riemann Zeta Function: Constraint Metric Tensor
             and the Bound $\Lambda \le 0.081$},
  author  = {Ismail, Mesut},
  year    = {2026},
  note    = {Submitted to J. Number Theory (JNTH-D-26-00613), v13}
}
```

## License

MIT

