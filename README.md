# RH Partial-Fraction Numerical Verifications

Code for reproducing numerical results in:

**"Partial-Fraction Constraints on Hypothetical Off-Line Zeros
of the Riemann Zeta Function"** -- Version 11, Mesut Ismail, March 2026

## Requirements

- Python 3.8+
- mpmath >= 1.3: `pip install mpmath`

## Usage

```bash
python run_all.py             # Full verification (~30 min with 200 zeros)
python run_all.py --quick     # Quick mode (~10 min)

# Individual verifications:
python verify_01_fundamental_identity.py
python verify_02_curvature_test.py
# ... etc.
```

## Architecture

Zeros are loaded in batches of 200 via `mpmath.zetazero(n)` at 25-digit
precision and cached to `zeros_cache.pkl`. The default loads 2000 zeros
(10 batches); adjust `N_BATCHES` in `config.py` for faster/slower runs.

## File Map

| File | Section | What it verifies |
|------|---------|-----------------|
| config.py | -- | Shared config, zero loading, Son computation |
| verify_01_* | 17.1 | V'(h,t) > 0 at 70 grid points (Thm 3.1) |
| verify_02_* | 17.2 | Curvature residual T_30, S/N > 10^4 (Thm 10.1) |
| verify_03_* | 17.3 | Sigma-equivalents >> 1 (Prop 11.1, Table 2) |
| verify_04_* | 17.4 | |Phi_off| = 1 on CL exactly (Thm 12.1) |
| verify_05_* | 17.5 | Concavity d^2/dt^2 < 0, 1 turning point/gap (Thm 14.1) |
| verify_06_* | 17.6 | Tunneling Mh0 = pi*s/sqrt(Cn), Cn stats (Thm 15.1) |
| verify_07_* | 17.7 | Self-consistency h0^2*Son = 2 to machine precision |
| verify_08_* | 17.8 | Speiser depth scaling, sqrt(delta) law (Thm 13.1) |
| verify_09_* | 17.9 | Collision time ODE, universal eta ~ 0.459 (Thm 16.1) |
| verify_10_* | 18-27 | Gram tensor E(T), budget B(T), E/B -> inf (Thm 17.2) |
| run_all.py | -- | Master runner |

## Verified Results (200 zeros, gamma <= 396)

| Verification | Result | Paper |
|-------------|--------|-------|
| V' > 0 | All 70 points | All 70 points |
| Curvature S/N | > 13,000 | > 10^4 |
| |Phi_off| on CL | 1.000000000000000 | 1 exact |
| Concavity | All 30 gaps, 1 TP each | All 1999 gaps |
| Cn range | [1.03, 1.76], mean 1.27 | [1.001, 2.306], mean 1.28 |
| h0^2 * Son | 2.0000000000 | 2.000000 +/- 1e-14 |
| Mh0/s | 2.80 +/- 0.18 | 2.80 +/- 0.20 |
| eta (collision) | 0.456 +/- 0.012 | 0.459 +/- 0.003 |
| E/B | grows as log^5(T) | Theorem 17.2 |
| Lambda (uncond.) | <= 0.200 | <= 0.200 |
| Lambda (cond.) | <= 0.047 | <= 0.047 |

## Notes

- With 200 zeros, Son ~ 0.22 for gap 1 vs paper's 3.09 with 2000 zeros.
  This does NOT affect structural results since h0 = sqrt(2/Son) adjusts
  accordingly and h0^2*Son = 2 holds exactly by construction.
- Collision time eta converges to ~0.456 with 200 zeros; with 2000 zeros
  and full ODE integration the paper obtains 0.459 +/- 0.003.
- The Speiser scaling constant C converges toward 1.044 with more zeros.
