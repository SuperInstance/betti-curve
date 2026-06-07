# betti-curve

> **Betti curves and Euler characteristic curves for persistent homology**

[![crates.io](https://img.shields.io/crates/v/betti-curve.svg)](https://crates.io/crates/betti-curve)
[![license](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

Betti curves track how the Betti numbers βₖ change as a function of the filtration parameter. Together with the Euler curve (alternating sum), they provide a compact summary of the topological evolution of a dataset.

## What It Computes

- **BettiCurve**: βₖ(ε) as ε increases — how many k-dimensional holes exist
- **BettiBarcode**: Interval representation of each topological feature's lifetime
- **BettiSummary**: Max, mean, total persistence statistics
- **EulerCurve**: χ(ε) = β₀ - β₁ + β₂ - ... at each scale
- **PersistenceEntropy**: Information-theoretic summary of barcode complexity

## Installation

```toml
[dependencies]
betti-curve = "0.1.0"
```

## License

MIT © [SuperInstance](https://github.com/SuperInstance)

---

*Part of the [Exocortex](https://github.com/SuperInstance/exocortex) project — persistent cognitive substrate for multi-agent systems.*
