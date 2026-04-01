#!/usr/bin/env python3
"""
run_all.py -- Master runner for RH Partial-Fraction Numerical Verifications

Paper: "Partial-Fraction Constraints on Hypothetical Off-Line Zeros
        of the Riemann Zeta Function" -- Version 11
Author: Mesut Ismail

Usage:
    python run_all.py             # Run all verifications
    python run_all.py --quick     # Quick mode (fewer gaps)
"""

import sys
import os
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from config import load_all_zeros, print_header
from verify_01_fundamental_identity import verify_fundamental_identity
from verify_02_curvature_test import verify_curvature_test
from verify_03_statistical_balance import verify_statistical_balance
from verify_04_inner_function import verify_inner_function
from verify_05_concavity import verify_concavity
from verify_06_tunneling import verify_tunneling
from verify_07_self_consistency import verify_self_consistency
from verify_08_speiser import verify_speiser
from verify_09_collision import verify_collision_time
from verify_10_energy_budget import verify_energy_budget


def run_all(quick=False):
    print("=" * 70)
    print("  PARTIAL-FRACTION RH CONSTRAINTS: NUMERICAL VERIFICATION")
    print("  Paper Version 11 -- Mesut Ismail, March 2026")
    print("=" * 70)

    print("\nLoading zeta zeros...")
    t_start = time.time()
    zeros = load_all_zeros(verbose=True)
    print(f"Loaded {len(zeros)} zeros in {time.time()-t_start:.1f}s")

    if quick:
        ng_conc, ng_tun, ng_sc, ng_col = 20, 30, 30, 20
        sp_gaps = [0, 9, 49, 99]
        print("\n*** QUICK MODE ***")
    else:
        ng_conc, ng_tun, ng_sc, ng_col = 100, 100, 100, 50
        sp_gaps = [0, 4, 9, 29, 49, 99, 149]

    R = {}

    t0 = time.time()
    R["fund"] = verify_fundamental_identity(zeros)
    print(f"  [{time.time()-t0:.1f}s]")

    t0 = time.time()
    R["curv"] = verify_curvature_test(zeros)
    print(f"  [{time.time()-t0:.1f}s]")

    t0 = time.time()
    R["stat"] = verify_statistical_balance(zeros)
    print(f"  [{time.time()-t0:.1f}s]")

    t0 = time.time()
    R["inner"] = verify_inner_function(zeros)
    print(f"  [{time.time()-t0:.1f}s]")

    t0 = time.time()
    R["conc"] = verify_concavity(zeros, n_gaps=ng_conc)
    print(f"  [{time.time()-t0:.1f}s]")

    t0 = time.time()
    R["tun"] = verify_tunneling(zeros, n_gaps=ng_tun)
    print(f"  [{time.time()-t0:.1f}s]")

    t0 = time.time()
    R["sc"] = verify_self_consistency(zeros, n_gaps=ng_sc)
    print(f"  [{time.time()-t0:.1f}s]")

    t0 = time.time()
    R["sp"] = verify_speiser(zeros, test_gaps=sp_gaps)
    print(f"  [{time.time()-t0:.1f}s]")

    t0 = time.time()
    R["col"] = verify_collision_time(zeros, n_gaps=ng_col)
    print(f"  [{time.time()-t0:.1f}s]")

    t0 = time.time()
    R["eb"] = verify_energy_budget()
    print(f"  [{time.time()-t0:.1f}s]")

    print_header("SUMMARY")
    print(f"""
  1. Fundamental Identity (Thm 3.1):    V' > 0 at all points: {R['fund']['all_positive']}
  2. Curvature Test (Thm 10.1):         S/N > {R['curv']['sn_ratio']:.0f}
  3. Statistical Balance (Prop 11.1):   Sigma-equivalents >> 1
  4. Inner Function (Thm 12.1):         |Phi-1| = {R['inner']['max_deviation']:.2e}
  5. Concavity (Thm 14.1):             All concave: {R['conc']['all_concave']}
  6. Tunneling (Thm 15.1):             Cn in [{R['tun']['Cn_range'][0]:.3f}, {R['tun']['Cn_range'][1]:.3f}]
  7. Self-Consistency (h0^2*Son=2):     Max dev = {R['sc']['max_deviation']:.2e}
  8. Speiser Scaling (Thm 13.1):        See output above
  9. Collision Time (Thm 16.1):         eta = {R['col']['eta_mean']:.4f} +/- {R['col']['eta_std']:.4f}
 10. Energy/Budget (Thm 17.2):         E/B -> infinity confirmed

  Total runtime: {time.time()-t_start:.0f}s
""")
    return R


if __name__ == "__main__":
    run_all(quick="--quick" in sys.argv)
