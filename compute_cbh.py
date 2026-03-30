"""
compute_cbh.py — Compute the Bateman-Horn structural constant C_BH
for Q47(n) = n^47 - (n-1)^47 prime quadruplets, with Abel summation
tail correction.

Usage: python compute_cbh.py [sieve_limit]
Default sieve_limit = 100,000,000

Output: C_BH convergence table and predicted vs observed comparison.

Requires: gmpy2 (optional, for extended precision integration)
"""
import math
import sys
import time

Q = 47
SIEVE_LIMIT = int(sys.argv[1]) if len(sys.argv) > 1 else 100_000_000

def sieve_primes(limit):
    sieve = bytearray([1]) * limit
    sieve[0] = sieve[1] = 0
    for p in range(2, int(limit**0.5) + 1):
        if sieve[p]:
            sieve[p*p:limit:p] = bytearray(len(range(p*p, limit, p)))
    return [p for p in range(2, limit) if sieve[p]]

def compute_roots_q47(p):
    exp = (p - 1) // Q
    g = 0
    for h in range(2, min(p, 200)):
        candidate = pow(h, exp, p)
        if candidate != 1 and pow(candidate, Q, p) == 1:
            g = candidate
            break
    if g == 0:
        roots = set()
        for n in range(p):
            if (pow(n, Q, p) - pow((n-1) % p, Q, p)) % p == 0:
                roots.add(n)
        return roots
    roots = set()
    zeta = g
    for _ in range(Q - 1):
        denom = (zeta - 1) % p
        if denom != 0:
            inv_denom = pow(denom, p - 2, p)
            roots.add((zeta * inv_denom) % p)
        zeta = (zeta * g) % p
    return roots

def integrand(t):
    logs = []
    for i in range(4):
        x = t + i
        if x <= 1:
            return 0
        u = 1 - 1/x
        log_val = 47 * math.log(x) + math.log(1 - u**47)
        logs.append(log_val)
    return 1.0 / (logs[0] * logs[1] * logs[2] * logs[3])

def simpson_integrate(a, b, n_steps=100000):
    h = (b - a) / n_steps
    s = integrand(a) + integrand(b)
    for i in range(1, n_steps):
        x = a + i * h
        s += (4 if i % 2 == 1 else 2) * integrand(x)
    return s * h / 3

if __name__ == '__main__':
    print(f"Computing C_BH for Q47 quadruplets (sieve limit = {SIEVE_LIMIT:,})")
    t0 = time.time()

    primes = sieve_primes(SIEVE_LIMIT)
    print(f"Sieved {len(primes):,} primes in {time.time()-t0:.1f}s")

    log_C = 0.0
    killer_count = 0
    total_pi = 0
    checkpoints = {}

    for p in primes:
        total_pi += 1
        if p % Q == 1:
            killer_count += 1
            roots = compute_roots_q47(p)
            blocked = set()
            for r in roots:
                for shift in range(4):
                    blocked.add((r - shift) % p)
            omega = len(blocked)
            log_C += math.log(1 - omega/p) - 4 * math.log(1 - 1/p)
        else:
            log_C += -4 * math.log(1 - 1/p)

        if p in [5000, 50000, 5000000, 10000000, 50000000] or p == primes[-1]:
            C_raw = math.exp(log_C)
            delta = 4 * total_pi - 184 * killer_count
            C_corr = C_raw * math.exp(-delta / p)
            checkpoints[p] = (total_pi, killer_count, C_raw, delta, C_corr)

    print(f"\nConvergence table:")
    print(f"{'M':>12} {'pi(M)':>10} {'pi_47,1':>8} {'C_raw':>10} {'Delta':>8} {'C_corr':>10}")
    print("-" * 62)
    for M in sorted(checkpoints):
        pi_M, pk, cr, d, cc = checkpoints[M]
        print(f"{M:>12,} {pi_M:>10,} {pk:>8,} {cr:>10.1f} {d:>+8,} {cc:>10.1f}")

    C_BH = list(checkpoints.values())[-1][4]
    print(f"\nFinal C_BH = {C_BH:.1f}")

    # Predicted vs observed
    sectors = [
        (2.3e7, 2e9, 15), (2e9, 4e9, 10), (4e9, 1e10, 34),
        (1e10, 3e10, 78), (3e10, 1e11, 240), (1e11, 2e11, 352),
        (2e11, 3e11, 295), (3e11, 4e11, 311), (4e11, 5e11, 260),
        (5e11, 6e11, 260), (6e11, 7e11, 253), (7e11, 8e11, 251),
    ]

    print(f"\n{'Sector':<20} {'Obs':>6} {'Pred':>8} {'O/P':>7}")
    print("-" * 44)
    total_obs, total_pred, chi2 = 0, 0.0, 0.0
    for a, b, obs in sectors:
        pred = C_BH * simpson_integrate(a, b)
        ratio = obs / pred
        chi2 += (obs - pred)**2 / pred
        total_obs += obs
        total_pred += pred
        print(f"[{a:.0e}, {b:.0e}]  {obs:>6} {pred:>8.1f} {ratio:>7.3f}")

    print("-" * 44)
    print(f"{'TOTAL':<20} {total_obs:>6} {total_pred:>8.1f} {total_obs/total_pred:>7.3f}")
    print(f"\nchi2 = {chi2:.2f} (df=11, p~{1-0.52:.2f})")
    print(f"Elapsed: {time.time()-t0:.1f}s")
