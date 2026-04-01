"""
Part II: Gram Tensor and Energy-Budget Ratio (Sections 18-27)

Two E/B bounds from Theorem 21.2:
  (i)  Speiser-free: E/B >= G_{II}/B = pi*log^3(T)/2  [valid at ALL depths]
  (ii) Full bound:   E/B >= c0*log^5(T), c0 = 120/(2pi)^5  [Theorem 17.2]

The log^5 bound is the headline result; it is dominated by G_{IV,IV} (Speiser class).
"""

import math
from config import print_header


def gram_diagonal(logT):
    """Gram tensor diagonal elements at height T (using logT directly)."""
    two_pi = 2 * math.pi
    return {
        "I":   logT**2 / (4 * math.pi**2),
        "II":  logT**4 / 4,                       # class II: b=1/2, h0-independent
        "III": 6 * logT**4 / two_pi**4,
        "IV":  120 * logT**6 / two_pi**6,          # dominant term -> log^5 in E/B
        "V":   math.log(logT / two_pi) if logT > two_pi else 0.01
    }


def verify_energy_budget(verbose=True):
    if verbose:
        print_header("Part II: Energy-Budget Ratio (Sections 18-27)")

    log10_values = [12, 23, 100, 1000]

    # --- Both bounds ---
    if verbose:
        print("\nTheorem 21.2: Two E/B bounds")
        print(f"\n  {'T':>10s}  {'G_II/B':>12s}  {'pi*log^3/2':>12s}  "
              f"{'E/B (full)':>12s}  {'c0*log^5':>12s}")
        print(f"  {'-'*10}  {'-'*12}  {'-'*12}  {'-'*12}  {'-'*12}")

    c0 = 120 / (2 * math.pi)**5
    ratios = []

    for log10T in log10_values:
        logT = log10T * math.log(10)
        G = gram_diagonal(logT)
        E = sum(G.values())
        B = logT / (2 * math.pi)

        ratio_II = G["II"] / B          # = pi*log^3/2
        sf_pred  = math.pi * logT**3 / 2

        ratio_full = E / B
        c0_pred    = c0 * logT**5

        ratios.append((log10T, ratio_full))

        if verbose:
            print(f"  10^{log10T:>4d}  {ratio_II:12.2e}  {sf_pred:12.2e}  "
                  f"{ratio_full:12.2e}  {c0_pred:12.2e}")

    if verbose:
        print(f"\n  c0 = 120/(2*pi)^5 = {c0:.6f}")
        print(f"\n  (i)  G_II/B  = pi*log^3(T)/2   -- valid at ALL depths incl. self-consistency")
        print(f"  (ii) E/B     >= c0*log^5(T)     -- Theorem 17.2 (outside SC regime)")
        print(f"  Both -> infinity as T -> infinity.")

    # --- Individual condition costs ---
    if verbose:
        print(f"\nIndividual condition costs at T = 10^12 (Section 25):")
        logT = 12 * math.log(10)
        G = gram_diagonal(logT)
        h0 = math.pi / logT
        log_h0 = abs(math.log(h0))
        B = logT / (2 * math.pi)

        for name, key, scaling in [
            ("(a) Apollonius", "I",   "log^2 T"),
            ("(b) Sign barrier", "II",  "log^4 T"),
            ("(c) SC / (e) Curv", "III", "log^4 T"),
            ("(h) Speiser",     "IV",  "log^6 T"),
            ("(g) Inner/Euler", "V",   "log log T"),
        ]:
            eps_k = G[key] / log_h0
            print(f"    {name:>20s} (class {key:>3s}):  eps = {eps_k:12.1f}   [{scaling}]")
        print(f"    {'Budget':>20s}:               B = {B:12.2f}   [log T]")
        print(f"\n    Speiser (h) alone exceeds budget by factor ~{G['IV']/log_h0/B:.0f}")

    # --- Lambda bounds ---
    if verbose:
        print(f"\nLambda bounds (Theorem 17.1):")
        Lmax = 1.614
        # (i) Unconditional
        eta2 = 0.615
        h0_unc = Lmax / 2  # Cn >= 1
        L_unc = eta2 * h0_unc**2 / 2
        print(f"  (i)   Unconditional: Lambda <= {eta2}*{h0_unc:.3f}^2/2 = {L_unc:.3f}")
        # (ii) Conditional
        h0_cond = Lmax / (2 * math.sqrt(6.95))
        L_cond = h0_cond**2 / 2
        print(f"  (ii)  Conditional (Cn>=6.95): h0 <= {h0_cond:.3f}, Lambda <= {L_cond:.3f}")
        # (iii) With eta
        eta = 0.459
        L_eta = eta * h0_cond**2 / 2
        print(f"  (iii) With eta=0.459: Lambda <= {L_eta:.3f}")

    return {"ratios": ratios, "c0": c0}


if __name__ == "__main__":
    verify_energy_budget()
