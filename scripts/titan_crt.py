"""
titan_crt.py — CRT-targeted search for 1000-digit Q47 prime quadruplets.

Uses the Chinese Remainder Theorem to construct starting coordinates
that bypass the strongest congruent primes, then steps through
candidates testing for quadruplets of probable primes.

Requires: gmpy2

Usage: python titan_crt.py [target_digits]
  Default: target_digits=1000
"""
import time, sys, random
from math import gcd

try:
    import gmpy2
except ImportError:
    print("Error: gmpy2 required. Install via: pip install gmpy2")
    sys.exit(1)

sys.set_int_max_str_digits(2_000_000)

Q = 47
TARGET_DIGITS = int(sys.argv[1]) if len(sys.argv) > 1 else 1000

def q47(n):
    return n**Q - (n-1)**Q

def is_prp(q_val):
    return gmpy2.powmod(2, q_val - 1, q_val) == 1

def get_congruent_primes(count):
    primes = []
    p = 2
    while len(primes) < count:
        if gmpy2.is_prime(p) and p % Q == 1:
            primes.append(p)
        p += 1
    return primes

def get_safe_residues(p):
    blocked = set()
    for x in range(p):
        if (pow(x, Q, p) - pow(x-1, Q, p)) % p == 0:
            for shift in range(4):
                blocked.add((x - shift) % p)
    return [x for x in range(p) if x not in blocked]

def extended_gcd(a, b):
    if a == 0:
        return b, 0, 1
    g, y, x = extended_gcd(b % a, a)
    return g, x - (b // a) * y, y

def crt(remainders, moduli):
    n_sum = 0
    M = 1
    for m in moduli:
        M *= m
    for r, m in zip(remainders, moduli):
        p = M // m
        _, inv, _ = extended_gcd(p, m)
        n_sum += r * inv * p
    return n_sum % M, M

if __name__ == '__main__':
    target_val = 10**(TARGET_DIGITS - 1) // 47
    base_n, _ = gmpy2.iroot(gmpy2.mpz(target_val), 46)
    base_n = int(base_n)
    print(f"Target: {TARGET_DIGITS}-digit primes")
    print(f"Base coordinate: {base_n} ({len(str(base_n))} digits)")

    K = 10
    print(f"Computing safe residues for {K} congruent primes...")
    killers = get_congruent_primes(K)
    safe = {p: get_safe_residues(p) for p in killers}
    chosen = [random.choice(safe[p]) for p in killers]
    n_crt, M = crt(chosen, killers)

    offset = (n_crt - base_n) % M
    n_start = base_n + offset
    print(f"CRT step size M = {M}")
    print(f"Searching...")

    attempts = 0
    t0 = time.time()
    current = n_start

    while True:
        attempts += 1
        p0 = q47(current)
        if is_prp(p0):
            p1 = q47(current + 1)
            if is_prp(p1):
                p2 = q47(current + 2)
                if is_prp(p2):
                    p3 = q47(current + 3)
                    if is_prp(p3):
                        ts = time.strftime('%Y-%m-%d %H:%M:%S')
                        print(f"\n[{ts}] FOUND QUADRUPLET at N = {current}")
                        print(f"  Q47(N) has {len(str(p0))} digits")
                        with open(f"Q47_crt_{TARGET_DIGITS}d.log", "a") as f:
                            f.write(f"[{ts}] N = {current}\n")
                            for i in range(4):
                                f.write(f"  ({current}+{i})^47-({current}+{i}-1)^47\n")
                            f.write("\n")

        current += M
        if attempts % 5000 == 0:
            rate = attempts / (time.time() - t0)
            sys.stdout.write(f"\r  Attempts: {attempts:,} | Rate: {rate:.1f}/s ")
            sys.stdout.flush()
