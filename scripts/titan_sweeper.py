"""
titan_sweeper.py — Modular sieve for Q47(n) prime quadruplets.

Searches a range [START_N, END_N] for four consecutive integers n, n+1, n+2, n+3
such that Q47(n+i) = (n+i)^47 - (n+i-1)^47 is a probable prime for all i in {0,1,2,3}.

Requires: gmpy2

Usage: python titan_sweeper.py [start] [end]
  Default: start=30_000_000_000, end=100_000_000_000
"""
import time, sys, multiprocessing, os

try:
    import gmpy2
except ImportError:
    print("Error: gmpy2 required. Install via: pip install gmpy2")
    sys.exit(1)

Q = 47
START_N = int(float(sys.argv[1])) if len(sys.argv) > 1 else 30_000_000_000
END_N   = int(float(sys.argv[2])) if len(sys.argv) > 2 else 100_000_000_000
CHUNK_SIZE = 10_000_000
SIEVE_LIMIT = 5_000_000

def q47(n):
    return n**Q - (n-1)**Q

def generate_primes(limit):
    sieve = [True] * limit
    for p in range(2, int(limit**0.5) + 1):
        if sieve[p]:
            for i in range(p*p, limit, p):
                sieve[i] = False
    return [p for p in range(2, limit) if sieve[p]]

def precompute_invalid_residues(primes):
    invalid_map = {}
    for p in primes:
        if p % Q != 1:
            continue
        roots = []
        for x in range(p):
            if pow(x, Q, p) == pow(x-1, Q, p):
                roots.append(x)
        if roots:
            blocked = set()
            for r in roots:
                for shift in range(4):
                    blocked.add((r - shift) % p)
            invalid_map[p] = blocked
    return invalid_map

def worker(start_n, chunk_size, invalid_map, result_queue):
    actual = min(chunk_size, END_N - start_n)
    if actual <= 0:
        return
    candidates = bytearray([1]) * actual
    for p, blocked in invalid_map.items():
        for inv_r in blocked:
            rem = start_n % p
            first = inv_r - rem if rem <= inv_r else inv_r - rem + p
            candidates[first::p] = b'\x00' * len(candidates[first::p])

    for i in range(actual - 3):
        if not candidates[i]:
            continue
        n = start_n + i
        if not gmpy2.is_prime(q47(n), 25): continue
        if not gmpy2.is_prime(q47(n+1), 25): continue
        if not gmpy2.is_prime(q47(n+2), 25): continue
        if not gmpy2.is_prime(q47(n+3), 25): continue
        result_queue.put(n)
    result_queue.put(('DONE', start_n, actual))

def main():
    print(f"Q47 Sweeper: [{START_N:,}, {END_N:,}]")
    primes = generate_primes(SIEVE_LIMIT)
    invalid_map = precompute_invalid_residues(primes)
    print(f"Sieve primes: {len(invalid_map)} congruent primes loaded")

    result_queue = multiprocessing.Queue()
    num_workers = os.cpu_count() or 4
    current = START_N
    active = []
    found = 0
    t0 = time.time()

    logfile = f"Q47_sweep_{START_N}_{END_N}.log"

    while current < END_N or active:
        while len(active) < num_workers and current < END_N:
            p = multiprocessing.Process(target=worker,
                args=(current, CHUNK_SIZE, invalid_map, result_queue))
            p.start()
            active.append(p)
            current += CHUNK_SIZE

        for p in active[:]:
            if not p.is_alive():
                p.join()
                active.remove(p)

        while not result_queue.empty():
            msg = result_queue.get()
            if isinstance(msg, tuple) and msg[0] == 'DONE':
                pct = min(100, 100*(msg[1]-START_N)/(END_N-START_N))
                sys.stdout.write(f"\rProgress: {pct:.1f}%  Found: {found}")
                sys.stdout.flush()
            else:
                n = msg
                found += 1
                ts = time.strftime('%Y-%m-%d %H:%M:%S')
                print(f"\n[{ts}] QUADRUPLET: {n}, {n+1}, {n+2}, {n+3}")
                with open(logfile, "a") as f:
                    f.write(f"[{ts}] FOUND QUADRUPLET: {n}, {n+1}, {n+2}, {n+3}\n")
        time.sleep(0.1)

    print(f"\nDone. Found {found} quadruplets in {(time.time()-t0)/3600:.2f} hours")

if __name__ == "__main__":
    multiprocessing.freeze_support()
    main()
