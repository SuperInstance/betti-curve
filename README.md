# betti-curve

> **Betti curves and Euler characteristic curves — track the topological evolution of your data across scales**

[![crates.io](https://img.shields.io/crates/v/betti-curve.svg)](https://crates.io/crates/betti-curve)
[![docs.rs](https://docs.rs/betti-curve/badge.svg)](https://docs.rs/betti-curve)
[![license](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

## What is a Betti Curve?

In persistent homology, topological features (connected components, loops, voids) are born and die as you vary a scale parameter ε. A **barcode** records each feature's lifetime as an interval [b, d). A **Betti curve** βₖ(ε) counts how many k-dimensional features are alive at each scale — it's a step function that rises when features are born and falls when they die.

The **Euler characteristic curve** χ(ε) = β₀(ε) − β₁(ε) + β₂(ε) − ... compresses all dimensions into a single integer-valued function that captures the overall topological complexity at each scale.

Together with **persistence entropy** (a Shannon entropy measure over bar lengths), these curves provide functional summaries of topological data suitable for statistical analysis, machine learning features, and visual interpretation.

## Why Does This Matter?

Betti curves and Euler curves are among the most useful summaries in topological data analysis:

- **Functional data analysis**: Betti curves are functions you can feed into FDA methods — smoothing, PCA, regression
- **Scale selection**: Peaks in Betti curves indicate "interesting" scales where topology is richest
- **Classification**: The shape of Betti curves distinguishes different data-generating processes
- **Complexity monitoring**: Euler curves track total topological complexity as a single number
- **Entropy**: Persistence entropy quantifies the diversity of topological feature lifetimes

Real-world applications:
- **Protein folding**: Track how secondary structure elements (loops, tunnels) appear during folding
- **Network analysis**: Monitor connected components and cycles in dynamic graphs
- **Image analysis**: Characterize texture via the topological signature across scales
- **Cosmology**: Study the topology of the cosmic web (voids, filaments, clusters)

## Architecture

```
┌───────────────────────────────────────────────────────────────┐
│                   Betti Curve Pipeline                         │
│                                                               │
│  Barcodes (H₀, H₁, H₂, ...)                                  │
│  ┌─────────────────────┐                                      │
│  │ H₀: ████  ██████    │    BettiCurve                        │
│  │ H₁:    ████         │──▶ β₀(ε): ───┐   ┌────             │
│  │ H₂:       ██        │    β₁(ε): ────┘   └───             │
│  └─────────────────────┘                                      │
│          │                                                    │
│          ▼                                                    │
│  ┌─────────────────┐  ┌──────────────────┐  ┌─────────────┐  │
│  │  EulerCurve     │  │  BettiSummary    │  │ Persistence │  │
│  │ χ(ε) = β₀-β₁+β₂│  │  max, mean,      │  │ Entropy     │  │
│  │                 │  │  total persist.  │  │ H = -Σp·ln p│  │
│  └─────────────────┘  └──────────────────┘  └─────────────┘  │
│                                                               │
│  ┌──────────────────────────────────────────────────────────┐ │
│  │  β(ε)  ▲                                                 │ │
│  │  3 │   ┤─────┐                                           │ │
│  │  2 │   │     └──────┐                                    │ │
│  │  1 │   │            └─────                               │ │
│  │  0 └───┴──────────────────┴─────▶ ε                      │ │
│  └──────────────────────────────────────────────────────────┘ │
└───────────────────────────────────────────────────────────────┘
```

## Quick Start

```rust
use betti_curve::{BettiBarcode, BettiCurve, BettiSummary, EulerCurve, PersistenceEntropy};

// Create a barcode for H₀ (connected components)
let barcode = BettiBarcode::new(0, vec![
    (0.0, 1.0),   // component 1 merges at ε=1.0
    (0.0, 3.0),   // component 2 merges at ε=3.0
    (0.5, 2.5),   // component 3 merges at ε=2.5
]);

// Build a Betti curve over [0, 3.5] with 100 samples
let curve = BettiCurve::from_barcode(&barcode, 0.0, 3.5, 100);
println!("Max Betti number: {}", curve.max_betti());
println!("Betti at ε=0.7: {}", curve.value_at(0.7));
println!("Area under curve: {:.2}", curve.area());

// Get summary statistics
let summary = BettiSummary::from_barcode(&barcode);
println!("Max: {}, Mean: {:.2}, Total persistence: {:.2}",
    summary.max, summary.mean, summary.total_persistence);
```

### Euler Characteristic Curve

```rust
// Create barcodes for multiple dimensions
let h0 = BettiBarcode::new(0, vec![(0.0, 5.0), (0.0, 5.0)]);
let h1 = BettiBarcode::new(1, vec![(1.0, 3.0)]);
let h2 = BettiBarcode::new(2, vec![]);

// Compute Euler curve: χ = β₀ − β₁ + β₂
let euler = EulerCurve::from_barcodes(&[h0, h1, h2], 0.0, 5.0, 50);
println!("χ range: [{}, {}]", euler.min(), euler.max());
```

### Persistence Entropy

```rust
let barcode = BettiBarcode::new(0, vec![
    (0.0, 2.0),
    (0.5, 4.0),
    (1.0, 1.5),
]);

// Entropy at a specific scale
let entropy_at_1 = PersistenceEntropy::entropy_at(&barcode, 1.0);

// Global entropy over all bars
let global = PersistenceEntropy::global_entropy(&barcode);
println!("Global persistence entropy: {:.4}", global);
```

## API Reference

### BettiBarcode

| Method | Returns | Description |
|--------|---------|-------------|
| `BettiBarcode::new(dim, bars)` | `BettiBarcode` | Create barcode for homology dimension `dim` |
| `barcode.alive_at(t)` | `usize` | Count features alive at parameter t |
| `barcode.total_persistence()` | `f64` | Sum of all bar lengths |
| `barcode.num_bars()` | `usize` | Number of intervals |
| `barcode.min_birth()` | `f64` | Earliest birth value |
| `barcode.max_death()` | `f64` | Latest death value |

### BettiCurve

| Method | Returns | Description |
|--------|---------|-------------|
| `BettiCurve::from_barcode(barcode, t_min, t_max, n)` | `BettiCurve` | Build curve with `n` samples |
| `curve.value_at(t)` | `usize` | Interpolated Betti number at t |
| `curve.max_betti()` | `usize` | Maximum Betti number observed |
| `curve.area()` | `f64` | Area under the curve (trapezoidal rule) |

### BettiSummary

| Method | Returns | Description |
|--------|---------|-------------|
| `BettiSummary::from_barcode(barcode)` | `BettiSummary` | Compute max, mean, total persistence |

### EulerCurve

| Method | Returns | Description |
|--------|---------|-------------|
| `EulerCurve::from_barcodes(barcodes, t_min, t_max, n)` | `EulerCurve` | χ(ε) from multiple dimensions |
| `euler.max()` | `i64` | Maximum Euler characteristic |
| `euler.min()` | `i64` | Minimum Euler characteristic |

### PersistenceEntropy

| Method | Returns | Description |
|--------|---------|-------------|
| `PersistenceEntropy::entropy_at(barcode, t)` | `f64` | Shannon entropy at scale t |
| `PersistenceEntropy::global_entropy(barcode)` | `f64` | Entropy over all bars |

## Mathematical Background

### Betti Numbers

The k-th Betti number βₖ counts the number of k-dimensional "holes" in a topological space:
- β₀ = number of connected components
- β₁ = number of loops (independent 1-cycles)
- β₂ = number of voids (cavities)

In persistent homology, Betti numbers vary with the filtration parameter ε:
```
βₖ(ε) = #{(bᵢ, dᵢ) : bᵢ ≤ ε < dᵢ}
```

### Euler Characteristic

The Euler characteristic is the alternating sum of Betti numbers:
```
χ(ε) = β₀(ε) − β₁(ε) + β₂(ε) − ...
```

This is a topological invariant that satisfies the Euler-Poincaré formula: for a simplicial complex with nₖ simplices of dimension k, χ = Σₖ(−1)ᵏnₖ = Σₖ(−1)ᵏβₖ.

### Persistence Entropy

Persistence entropy treats bar lengths as a probability distribution:
```
pᵢ = (dᵢ − bᵢ) / Σⱼ(dⱼ − bⱼ)
H = −Σᵢ pᵢ · ln(pᵢ)
```

- H = 0: all features have equal persistence (one bar, or all identical)
- H = ln(n): n features with equal persistence (maximum diversity)
- Higher entropy → more diverse feature lifetimes

### Area Under Betti Curve

The area Aₖ = ∫ βₖ(ε) dε is related to the total persistence: it measures the cumulative topological complexity across all scales. Larger area = more topological structure overall.

## Installation

```bash
cargo add betti-curve
```

Or add to your `Cargo.toml`:

```toml
[dependencies]
betti-curve = "0.1.0"
```

## Related Crates

- [`persistence-landscape`](https://github.com/SuperInstance/persistence-landscape) — Persistence landscapes for statistical TDA
- [`mapper-graph`](https://github.com/SuperInstance/mapper-graph) — Mapper algorithm for point cloud topology
- [`cech-complex`](https://github.com/SuperInstance/cech-complex) — Čech complex construction
- [`witness-complex`](https://github.com/SuperInstance/witness-complex) — Witness complex approximation

## License

MIT © [SuperInstance](https://github.com/SuperInstance)

---

*Part of the [Exocortex](https://github.com/SuperInstance/exocortex) project — persistent cognitive substrate for multi-agent systems.*
