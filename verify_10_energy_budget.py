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

    # --- Lambda bound (Theorem 17.1 / Corollary, v14) ---
    if verbose:
        print(f"\nLambda bound (Theorem 17.1):")
        Lmax = 1.614
        half_sq = (Lmax / 2.0)**2     # (L_n/2)^2 = 0.651...

        # The UNCONDITIONAL bound uses the critical-strip depth ceiling h0 < 1/2
        # (NOT the self-consistency depth 0.807), together with the rigorous
        # two-zero collision integral.  The full backward-flow collision time is
        # bounded above by this integral via the comparison principle.
        #   tau_2(h0=1/2, L_n=Lmax) = int_0^{1/2} h(h^2 + (L/2)^2)/((L/2)^2 + 5h^2) dh
        n = 200000
        a, b = 0.0, 0.5
        dx = (b - a) / n
        def integrand(h):
            return h * (h**2 + half_sq) / (half_sq + 5.0 * h**2)
        total = integrand(a) + integrand(b)
        for i in range(1, n):
            x = a + i * dx
            total += (4 if i % 2 else 2) * integrand(x)
        tau2 = total * dx / 3.0
        eta2_at_Cn1 = 0.487  # ratio tau2/(h0^2/2) at the self-consistent depth h0=L_n/2 (C_n=1)

        print(f"  Critical-strip depth ceiling (unconditional): h0 < 1/2")
        print(f"  Max gap at T0=3e12 (Trudgian): L_max = {Lmax}")
        print(f"  Rigorous two-zero integral at h0=1/2, L_n=L_max:")
        print(f"    Lambda <= tau_2(1/2, {Lmax}) = {tau2:.4f}  ->  Lambda <= 0.081")
        print(f"    (Lambda is the integral value at (h0,L_n)=(1/2,L_max), Layer A, rigorous;")
        print(f"     the dimensionless ratio eta_2 = {eta2_at_Cn1} refers to C_n=1 and is not used here.)")
        print(f"  This does NOT use the self-consistency depth (0.807) and does NOT")
        print(f"  use the numerical full-ODE eta ~ 0.46 (Layer B, approximately stable).")

    return {"ratios": ratios, "c0": c0}


if __name__ == "__main__":
    verify_energy_budget()
