"""
Section 17.8: Speiser Depth Scaling (Theorem 13.1)

At self-consistency (h0 = hthr), no Speiser zero exists (when using effective
threshold h_eff_thr based on Seff). For h0 > h_eff_thr:
  h' = C * sqrt(delta * h_eff_thr), C ~ 1.044 (paper).

The Speiser zero is where Re[xi'/xi + hypothetical_pole] = 0.
The hypothetical pole: 1/(h'-h0) + 1/(h'+h0) = 2h'/(h'^2 - h0^2).
"""

import mpmath
import math
import statistics
from config import print_header


def xi_prime_over_xi_re(h_val, t0):
    s = mpmath.mpc(0.5 + h_val, t0)
    psi = mpmath.digamma(s / 2)
    zv = mpmath.zeta(s)
    if abs(zv) < 1e-30:
        return None
    zp = mpmath.diff(mpmath.zeta, s)
    return float(mpmath.re(1/s + 1/(s-1) - mpmath.log(mpmath.pi)/2 + psi/2 + zp/zv))


def G_total(hp, t0, h0):
    """G(h') = Re[xi'/xi] + pole from hypothetical off-line zero pair."""
    base = xi_prime_over_xi_re(hp, t0)
    if base is None:
        return None
    if abs(hp**2 - h0**2) < 1e-20:
        return None
    pole = 2 * hp / (hp**2 - h0**2)
    return base + pole


def compute_Seff(t0, eps=1e-4):
    """Effective Speiser derivative via Richardson extrapolation."""
    g1 = xi_prime_over_xi_re(eps, t0)
    g2 = xi_prime_over_xi_re(2*eps, t0)
    d1 = g1 / eps
    d2 = g2 / (2*eps)
    return (4*d1 - d2) / 3


def find_speiser_zero(t0, h0, h_min=0.001):
    """Find h' where G_total(h', t0, h0) = 0 via bisection."""
    h_max = h0 - 0.001
    if h_max <= h_min:
        return None
    n = 300
    prev_g = None
    prev_h = None
    for j in range(n + 1):
        hp = h_min + (h_max - h_min) * j / n
        g = G_total(hp, t0, h0)
        if g is None:
            continue
        if prev_g is not None and prev_g * g < 0:
            a, b = prev_h, hp
            ga = prev_g
            for _ in range(60):
                mid = (a + b) / 2
                gm = G_total(mid, t0, h0)
                if gm is None:
                    break
                if ga * gm < 0:
                    b = mid
                else:
                    a = mid
                    ga = gm
            return (a + b) / 2
        prev_g = g
        prev_h = hp
    return None


def verify_speiser(zeros, test_gaps=None, verbose=True):
    if verbose:
        print_header("17.8 Speiser Depth Scaling (Theorem 13.1)")

    if test_gaps is None:
        test_gaps = [0, 4, 9, 29, 49, 99, 149]
    test_gaps = [g for g in test_gaps if g < len(zeros) - 1]

    if verbose:
        print("\nPart A: Effective threshold and Speiser test at self-consistency")
        print(f"  {'Gap':>5s}  {'t0':>8s}  {'Son':>8s}  {'Seff':>8s}  {'hthr':>8s}  {'h_eff':>8s}  {'SP@SC?':>8s}")

    results = []
    for idx in test_gaps:
        t0 = (float(zeros[idx]) + float(zeros[idx + 1])) / 2
        son = sum(1.0 / (t0 - float(gk))**2 for gk in zeros if abs(t0 - float(gk)) > 1e-12)
        seff = compute_Seff(t0)
        hthr = math.sqrt(2.0 / son)
        h_eff = math.sqrt(2.0 / seff) if seff > 0 else float('inf')

        sp_at_eff = find_speiser_zero(t0, h_eff)

        if verbose:
            print(f"  {idx+1:5d}  {t0:8.1f}  {son:8.4f}  {seff:8.4f}  {hthr:8.4f}  {h_eff:8.4f}  "
                  f"{'Yes' if sp_at_eff else 'No':>8s}")

        # Square-root scaling
        gap_C = []
        for delta in [0.05, 0.1, 0.2, 0.4]:
            h0 = h_eff + delta
            hp = find_speiser_zero(t0, h0)
            if hp:
                C = hp / math.sqrt(delta * h_eff)
                gap_C.append((delta, hp, C))

        results.append({
            "gap": idx + 1, "t0": t0, "son": son, "seff": seff,
            "hthr": hthr, "h_eff": h_eff,
            "sp_at_sc": sp_at_eff is not None,
            "scaling": gap_C
        })

    if verbose:
        print(f"\nPart B: Square-root scaling for h0 > h_eff_thr")
        print(f"  {'Gap':>5s}  {'delta':>8s}  {'h_prime':>10s}  {'C':>8s}")
        C_all = []
        for r in results:
            for delta, hp, C in r["scaling"]:
                C_all.append(C)
                print(f"  {r['gap']:5d}  {delta:8.3f}  {hp:10.6f}  {C:8.3f}")

        if C_all:
            print(f"\n  C = {statistics.mean(C_all):.3f} +/- {statistics.stdev(C_all):.3f}")

        print(f"\n  Paper reference: C = 1.044 +/- 0.001 (with 2000 zeros)")
        print(f"  Note: exact C depends on zero count; convergence improves with more zeros.")

    return results


if __name__ == "__main__":
    from config import load_all_zeros
    zeros = load_all_zeros()
    verify_speiser(zeros)
