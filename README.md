# Q47-Bateman-Horn

**Empirical Verification of the Bateman–Horn Conjecture for Q₄₇(n) = n⁴⁷ − (n−1)⁴⁷**

This repository contains the complete data, source code, ECPP certificates, and paper for a large-scale computational verification of the [Bateman–Horn conjecture](https://en.wikipedia.org/wiki/Bateman%E2%80%93Horn_conjecture) applied to the degree-46 polynomial Q₄₇(n) = n⁴⁷ − (n−1)⁴⁷ in the prime quadruplet configuration.

## Key Results

| Metric | Value |
|--------|-------|
| Prime quadruplets discovered | **2,359** |
| Search range (systematic sweep) | n ∈ [2.3×10⁷, 8×10¹¹] |
| ECPP-certified 1000-digit quadruplets | **11** (44 proven primes) |
| Bateman–Horn constant C_BH | **6514.2** (Abel-summed tail correction) |
| χ² goodness-of-fit | **10.56** (df=11, p ≈ 0.48) |
| Observed / Predicted ratio | **0.981** |

## Repository Structure

```
Q47-Bateman-Horn/
├── paper/
│   ├── Titanic_Quadruplet_Q47.tex    # LaTeX source
│   └── Titanic_Quadruplet_Q47.pdf    # Compiled paper
├── data/
│   ├── quadruplets_all_2359.csv      # All 2,359 starting coordinates
│   ├── ecpp_coordinates.csv          # 11 ECPP-certified coordinates
│   ├── bateman_horn_comparison.csv   # Observed vs predicted (12 sectors)
│   └── cbh_convergence.csv           # C_BH convergence table
├── scripts/
│   ├── compute_cbh.py                # Bateman-Horn constant computation
│   ├── titan_sweeper.py              # Modular sieve (systematic search)
│   └── titan_crt.py                  # CRT targeting (1000-digit search)
├── figures/
│   ├── figure1_OP_ratio.tex          # Figure 1 source (pgfplots)
│   └── figure1_OP_ratio.pdf          # Figure 1 compiled
├── ecpp_certificates/                # Primo ECPP certificate files
│   └── *.out                         # 48 certificate files (44 unique primes)
└── README.md
```

## The Polynomial

The polynomial Q₄₇(n) = n⁴⁷ − (n−1)⁴⁷ has degree 46 and possesses a remarkable **algebraic shield**: for any prime p ≢ 1 (mod 47), we have ω(p) = 0, meaning ~97.9% of all primes cannot divide any value of Q₄₇. This results in an unusually large Bateman–Horn structural constant (C_BH ≈ 6514, compared to ~1.32 for twin primes).

A **prime quadruplet** is a set of four consecutive integers {n, n+1, n+2, n+3} such that Q₄₇(n+i) is prime for all i ∈ {0,1,2,3}.

## Requirements

```bash
pip install gmpy2
```

For ECPP verification: [Primo 4.4.0](http://www.ellipsa.eu/public/primo/primo.html)

## Usage

### Systematic sweep (modular sieve)
```bash
# Search n ∈ [30 billion, 100 billion]
python scripts/titan_sweeper.py 3e10 1e11
```

### CRT targeting for 1000-digit primes
```bash
python scripts/titan_crt.py 1000
```

### Compute Bateman–Horn constant
```bash
python scripts/compute_cbh.py 100000000
```

## ECPP Certificates

All 44 constituent 1000-digit primes (from 11 quadruplets) have been rigorously proven prime using the Elliptic Curve Primality Proving (ECPP) algorithm as implemented in Primo 4.4.0. The full Atkin–Goldwasser–Kilian certificates are in `ecpp_certificates/`.

Each certificate file contains:
- `Status=Candidate certified prime` (the verdict)
- The full certificate chain reducing the 1000-digit input to a small prime

## Citation

```bibtex
@article{chen2026q47,
  author  = {Chen, Ruqing},
  title   = {Empirical Verification of the {Bateman--Horn} Conjecture
             for the Degree-46 Polynomial {$Q_{47}(n) = n^{47}-(n-1)^{47}$}},
  year    = {2026},
  note    = {\url{https://github.com/Ruqing1963/Q47-Bateman-Horn}}
}
```

## License

This work is released under the [MIT License](LICENSE).
