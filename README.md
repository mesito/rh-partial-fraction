# RH Partial-Fraction Numerical Verifications

Code for reproducing numerical results in:

**"Arithmetic Non-Realizability of Off-Line Blaschke Factors
in the Riemann Zeta Function: Constraint Metric Tensor and the Bound Λ ≤ 0.200"**
— Version 12, Mesut Ismail, May 2026

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
|------|---------|-----------------|
| `config.py` | — | Shared config, zero loading, S_on computation |
| `verify_01_*` | 17.1 | V'(h,t) > 0 at 70 grid points (Theorem 3.1) |
| `verify_02_*` | 17.2 | Curvature residual T₃₀, S/N > 10⁴ (Theorem 10.1) |
| `verify_03_*` | 17.3 | Sigma-equivalents ≫ 1 (Proposition 11.1) |
| `verify_04_*` | 17.4 | \|Φ_off\| = 1 on CL exactly (Theorem 12.1) |
| `verify_05_*` | 17.5 | Concavity d²/dt² < 0, 1 turning point/gap (Theorem 14.1) |
| `verify_06_*` | 17.6 | Tunneling Mh₀ = πs/√Cₙ, Cₙ statistics (Theorem 15.1) |
| `verify_07_*` | 17.7 | Self-consistency h₀²·S_on = 2 to machine precision |
| `verify_08_*` | 17.8 | Speiser depth scaling, √δ law (Theorem 13.1) |
| `verify_09_*` | 17.9 | Collision-time ODE, universal η ≈ 0.456 (Theorem 16.1) |
| `verify_10_*` | 18–27 | Gram tensor E(T), budget B(T), E/B → ∞ (Theorem 21.2) |
| `run_all.py` | — | Master runner |

## Verified Results (2000 zeros, γ ≤ 2515)

| Verification | Result (2000 zeros) | Paper claim |
|-------------|---------------------|-------------|
| V' > 0 | All 70 points ✓ | All 70 points |
| Curvature S/N | > 13,000 ✓ | > 10⁴ |
| \|Φ_off\| on CL | 1.000000000000000 ✓ | 1 exact |
| Concavity | All 50 gaps concave, 1 TP each ✓ | All gaps |
| Cₙ range | [1.002, 2.306], mean 1.29 ✓ | [1.002, 2.306], mean 1.29 |
| h₀² · S_on | 2.0000000000 ± 0 ✓ | 2.000000 ± 10⁻¹⁴ |
| Mh₀/s | 2.80 ± 0.19 ✓ | 2.80 ± 0.20 |
| Speiser C̄ | 1.18 ± 0.05 | 1.18 ± 0.05 |
| η (collision) | 0.456 ± 0.002 | 0.456 ± 0.002 |
| E/B | grows as log⁵ T ✓ | Theorem 21.2 |
| **Λ (unconditional)** | **≤ 0.200** ✓ | **≤ 0.200** |

## Key results

### Λ ≤ 0.200 (unconditional)

The main quantitative result. Combines three independent ingredients:

1. **Platt–Trudgian verification**: all zeros with γ ≤ 3 × 10¹² lie on the critical line.
2. **Self-consistency upper bound** (Theorem 5.3): h₀ ≤ Lₙ/(2√Cₙ) with Cₙ ≥ 1.
3. **Two-zero collision ODE** (Lemma in Section 16): η₂ = 0.615.

The two-zero case is a rigorous upper bound on collision time by comparison principle:
the full ODE (with all on-line zeros) has greater acceleration, hence faster collision.
This gives Λ ≤ η₂ · h₀²_max / 2 ≤ 0.615 × 0.807² / 2 ≈ 0.200.

This matches the Platt–Trudgian bound via an independent method — the first
new approach to the de Bruijn–Newman constant since Polymath15 (2019).

### Speiser depth scaling (Theorem 13.1)

At the self-consistent depth h₀ = √(2/S_on), no Speiser companion exists.
For h₀ > h_thr, a companion appears at depth h' ≈ C̄√(δ · h_thr) with
C̄ ≈ 1.18 ± 0.05 (square-root scaling, exponent α ≈ 0.51).

Note: gaps with large Lₙ (e.g., gap 1 with L₁ ≈ 6.89) have h_thr > ½,
so the Speiser companion cannot exist within the critical strip regardless
of h₀. The scaling law is verified on gaps with h_thr < ½ (e.g., gap 600).

### Structural inequality E/B → ∞ (Theorem 21.2)

The Speiser-free bound E/B ≥ π log³ T / 2 holds at all depths including
self-consistency, and is independent of the Speiser scaling constant.

## Errata from v11 → v12

| Item | v11 claim | v12 (corrected) | Impact |
|------|-----------|-----------------|--------|
| Selberg attribution | N₀(T) ≥ cT log T (Selberg 1942) | N₀(T) ≥ cT log T (**Levinson 1974**); Selberg proved N₀(T) ≥ cT | Attribution only; estimates unchanged |
| Lemma 13.2 (Caster drift) | S_drift ≈ −0.65 | **Removed** — drift is +0.012 (truncation residual converging to 0) | Caster framework removed |
| Speiser constant C | 1.044 ± 0.001 | **1.18 ± 0.05** | Quantitative only; phenomenology unchanged |
| Remark 5.4 (Cₙ range) | Cₙ ≥ 10.3 | **Cₙ ∈ [1.002, 2.306]**, mean 1.29 | Cn_min = 1.002, not ≥ 6.95 |
| Λ ≤ 0.047 conditional | Conditional on Cₙ ≥ 6.95 | **Removed** — condition not satisfied | Only Λ ≤ 0.200 unconditional remains |
| Λ ≤ 0.021 conditional | Conditional on η numerical | **Removed** | Only Λ ≤ 0.200 unconditional remains |
| η collision-time ratio | η = 0.459 ± 0.003 | **η = 0.456 ± 0.002** (300 nearby zeros per gap) | Numerical observation; Λ bound uses η₂ = 0.615 (analytical) |
| Speiser Step 4 verification | Gap 1 (S_on = 3.09, h_thr = 0.805) | **Gap 600** (S_on = 20.998, h_thr = 0.309) | Gap 1 has h_thr > ½, cannot test Speiser in critical strip |
| Gap-by-gap η table | S_on values from Lemma 13.2 | S_on from full 2000-zero computation | η ≈ 0.456 unchanged |

## Three-layer structure

The paper is organized into three layers with explicit logical status:

- **Layer A (rigorous):** 13 theorems proved from standard analytic number theory
  (Hadamard factorization, Levinson bound, Titchmarsh bounds, Ingham density estimate).
  Includes Λ ≤ 0.200 unconditional.

- **Layer B (numerical observations):** Gram tensor, collision-time universality,
  Speiser scaling. Structural inequality E/B → ∞ is rigorous; numerical constants
  (η = 0.456, C̄ = 1.18) are observations supported by 25-digit computation.

- **Layer C (conjectural):** Arithmetic Non-Realizability principle. Reformulates the
  obstruction at the prime side via Euler realizability, motivated by the Blaschke
  inner-function structure (Theorem 12.1).

## Citation

```bibtex
@article{ismail2026partialfraction,
  title   = {Arithmetic Non-Realizability of Off-Line Blaschke Factors
             in the Riemann Zeta Function: Constraint Metric Tensor
             and the Bound $\Lambda \le 0.200$},
  author  = {Ismail, Mesut},
  year    = {2026},
  note    = {Preprint, v12}
}
```

## License

MIT
