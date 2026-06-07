//! # betti-curve
//!
//! **Betti curves**, barcodes, Euler curves, and persistence entropy
//! for topological data analysis summaries.
//!
//! A Betti curve β_k(t) counts the number of k-dimensional topological
//! features alive at filtration parameter t. These curves provide a
//! functional summary of persistent homology suitable for statistical analysis.
//!
//! # Example
//!
//! ```
//! use betti_curve::{BettiBarcode, BettiCurve, BettiSummary};
//!
//! let barcode = BettiBarcode::new(0, vec![(0.0, 1.0), (0.5, 3.0), (2.0, 2.5)]);
//! let curve = BettiCurve::from_barcode(&barcode, 0.0, 3.0, 100);
//! assert!(curve.value_at(0.7) > 0);
//! let summary = BettiSummary::from_barcode(&barcode);
//! assert!(summary.total_persistence > 0.0);
//! ```

// ---------------------------------------------------------------------------
// BettiBarcode
// ---------------------------------------------------------------------------

/// A barcode for a single homology dimension: a list of (birth, death) intervals.
#[derive(Debug, Clone)]
pub struct BettiBarcode {
    /// Homology dimension.
    pub dim: usize,
    /// List of (birth, death) intervals.
    pub bars: Vec<(f64, f64)>,
}

impl BettiBarcode {
    /// Create a new barcode for homology dimension `dim`.
    pub fn new(dim: usize, bars: Vec<(f64, f64)>) -> Self {
        Self { dim, bars }
    }

    /// Number of bars.
    pub fn num_bars(&self) -> usize {
        self.bars.len()
    }

    /// True if the barcode is empty.
    pub fn is_empty(&self) -> bool {
        self.bars.is_empty()
    }

    /// Count bars alive at parameter `t`.
    pub fn alive_at(&self, t: f64) -> usize {
        self.bars
            .iter()
            .filter(|(b, d)| t >= *b && t < *d)
            .count()
    }

    /// Total persistence (sum of bar lengths).
    pub fn total_persistence(&self) -> f64 {
        self.bars.iter().map(|(b, d)| (d - b).max(0.0)).sum()
    }

    /// Maximum death value.
    pub fn max_death(&self) -> f64 {
        self.bars
            .iter()
            .map(|(_, d)| *d)
            .fold(f64::NEG_INFINITY, f64::max)
    }

    /// Minimum birth value.
    pub fn min_birth(&self) -> f64 {
        self.bars
            .iter()
            .map(|(b, _)| *b)
            .fold(f64::INFINITY, f64::min)
    }
}

// ---------------------------------------------------------------------------
// BettiCurve
// ---------------------------------------------------------------------------

/// A Betti curve: β(t) sampled at discrete points over a parameter range.
#[derive(Debug, Clone)]
pub struct BettiCurve {
    /// Parameter values (sorted).
    pub t_values: Vec<f64>,
    /// β(t) at each parameter value.
    pub betti_values: Vec<usize>,
    /// Homology dimension this curve represents.
    pub dim: usize,
}

impl BettiCurve {
    /// Build a Betti curve from a barcode over [t_min, t_max] with `num_samples` steps.
    pub fn from_barcode(barcode: &BettiBarcode, t_min: f64, t_max: f64, num_samples: usize) -> Self {
        let t_values: Vec<f64> = (0..=num_samples)
            .map(|i| t_min + (t_max - t_min) * i as f64 / num_samples as f64)
            .collect();
        let betti_values: Vec<usize> = t_values.iter().map(|&t| barcode.alive_at(t)).collect();
        Self {
            t_values,
            betti_values,
            dim: barcode.dim,
        }
    }

    /// Interpolate the Betti number at parameter t (step function: use nearest sample).
    pub fn value_at(&self, t: f64) -> usize {
        if self.t_values.is_empty() {
            return 0;
        }
        // Find nearest index
        let idx = self
            .t_values
            .binary_search_by(|probe| {
                probe
                    .partial_cmp(&t)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .unwrap_or_else(|i| i.min(self.betti_values.len() - 1));
        self.betti_values[idx]
    }

    /// Maximum Betti number in the curve.
    pub fn max_betti(&self) -> usize {
        self.betti_values.iter().copied().max().unwrap_or(0)
    }

    /// Area under the Betti curve (trapezoidal rule, as f64).
    pub fn area(&self) -> f64 {
        if self.t_values.len() < 2 {
            return 0.0;
        }
        let mut area = 0.0;
        for i in 0..self.t_values.len() - 1 {
            let dx = self.t_values[i + 1] - self.t_values[i];
            let avg = (self.betti_values[i] + self.betti_values[i + 1]) as f64 / 2.0;
            area += dx * avg;
        }
        area
    }
}

// ---------------------------------------------------------------------------
// BettiSummary
// ---------------------------------------------------------------------------

/// Statistical summary of a barcode.
#[derive(Debug, Clone)]
pub struct BettiSummary {
    /// Maximum Betti number observed.
    pub max: usize,
    /// Mean Betti number (averaged over sampled parameter range).
    pub mean: f64,
    /// Total persistence.
    pub total_persistence: f64,
}

impl BettiSummary {
    /// Compute summary statistics from a barcode.
    pub fn from_barcode(barcode: &BettiBarcode) -> Self {
        if barcode.is_empty() {
            return Self {
                max: 0,
                mean: 0.0,
                total_persistence: 0.0,
            };
        }
        let lo = barcode.min_birth();
        let hi = barcode.max_death();
        if hi <= lo {
            return Self {
                max: barcode.num_bars(),
                mean: barcode.num_bars() as f64,
                total_persistence: barcode.total_persistence(),
            };
        }
        let curve = BettiCurve::from_barcode(barcode, lo, hi, 200);
        Self {
            max: curve.max_betti(),
            mean: curve.betti_values.iter().sum::<usize>() as f64 / curve.betti_values.len() as f64,
            total_persistence: barcode.total_persistence(),
        }
    }
}

// ---------------------------------------------------------------------------
// EulerCurve
// ---------------------------------------------------------------------------

/// Euler characteristic curve: χ(t) = β₀(t) − β₁(t) + β₂(t) − …
#[derive(Debug, Clone)]
pub struct EulerCurve {
    pub t_values: Vec<f64>,
    pub euler_values: Vec<i64>,
}

impl EulerCurve {
    /// Compute an Euler curve from multiple barcodes of different dimensions.
    pub fn from_barcodes(barcodes: &[BettiBarcode], t_min: f64, t_max: f64, num_samples: usize) -> Self {
        let t_values: Vec<f64> = (0..=num_samples)
            .map(|i| t_min + (t_max - t_min) * i as f64 / num_samples as f64)
            .collect();
        let mut euler_values = vec![0i64; t_values.len()];

        for barcode in barcodes {
            let sign = if barcode.dim % 2 == 0 { 1i64 } else { -1i64 };
            for (i, &t) in t_values.iter().enumerate() {
                euler_values[i] += sign * barcode.alive_at(t) as i64;
            }
        }

        Self {
            t_values,
            euler_values,
        }
    }

    /// Value at the i-th sample.
    pub fn value_at_index(&self, i: usize) -> i64 {
        self.euler_values[i]
    }

    /// Maximum Euler characteristic.
    pub fn max(&self) -> i64 {
        *self.euler_values.iter().max().unwrap_or(&0)
    }

    /// Minimum Euler characteristic.
    pub fn min(&self) -> i64 {
        *self.euler_values.iter().min().unwrap_or(&0)
    }
}

// ---------------------------------------------------------------------------
// PersistenceEntropy
// ---------------------------------------------------------------------------

/// Information-theoretic summary based on persistence bar lengths.
pub struct PersistenceEntropy;

impl PersistenceEntropy {
    /// Compute persistence entropy at parameter t from a barcode.
    ///
    /// This treats the persistence of each bar alive at t as a probability
    /// distribution and computes the Shannon entropy.
    pub fn entropy_at(barcode: &BettiBarcode, t: f64) -> f64 {
        let persistences: Vec<f64> = barcode
            .bars
            .iter()
            .filter(|(b, d)| t >= *b && t < *d)
            .map(|(b, d)| (d - b).max(1e-15))
            .collect();
        if persistences.is_empty() {
            return 0.0;
        }
        let total: f64 = persistences.iter().sum();
        if total <= 0.0 {
            return 0.0;
        }
        persistences
            .iter()
            .map(|&p| {
                let prob = p / total;
                if prob > 0.0 { -prob * prob.ln() } else { 0.0 }
            })
            .sum()
    }

    /// Compute a single scalar entropy from the full barcode (ignoring time).
    pub fn global_entropy(barcode: &BettiBarcode) -> f64 {
        let persistences: Vec<f64> = barcode
            .bars
            .iter()
            .map(|(b, d)| (d - b).max(1e-15))
            .collect();
        if persistences.is_empty() {
            return 0.0;
        }
        let total: f64 = persistences.iter().sum();
        if total <= 0.0 {
            return 0.0;
        }
        persistences
            .iter()
            .map(|&p| {
                let prob = p / total;
                if prob > 0.0 { -prob * prob.ln() } else { 0.0 }
            })
            .sum()
    }
}

// ===========================================================================
// Tests
// ===========================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_barcode_alive_at() {
        let bc = BettiBarcode::new(0, vec![(0.0, 2.0), (1.0, 3.0)]);
        assert_eq!(bc.alive_at(0.5), 1);
        assert_eq!(bc.alive_at(1.5), 2);
        assert_eq!(bc.alive_at(2.5), 1);
        assert_eq!(bc.alive_at(3.5), 0);
    }

    #[test]
    fn test_barcode_total_persistence() {
        let bc = BettiBarcode::new(0, vec![(0.0, 2.0), (1.0, 4.0)]);
        assert!((bc.total_persistence() - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_barcode_min_max() {
        let bc = BettiBarcode::new(0, vec![(1.0, 5.0), (2.0, 8.0)]);
        assert!((bc.min_birth() - 1.0).abs() < 1e-10);
        assert!((bc.max_death() - 8.0).abs() < 1e-10);
    }

    #[test]
    fn test_betti_curve_from_barcode() {
        let bc = BettiBarcode::new(0, vec![(0.0, 1.0), (0.0, 2.0)]);
        let curve = BettiCurve::from_barcode(&bc, 0.0, 2.0, 4);
        assert_eq!(curve.dim, 0);
        assert_eq!(curve.t_values.len(), 5);
        assert!(curve.max_betti() >= 2);
    }

    #[test]
    fn test_betti_curve_area() {
        let bc = BettiBarcode::new(0, vec![(0.0, 10.0)]);
        let curve = BettiCurve::from_barcode(&bc, 0.0, 10.0, 100);
        // Area should be approximately 10.0
        assert!((curve.area() - 10.0).abs() < 0.5);
    }

    #[test]
    fn test_betti_curve_value_at() {
        let bc = BettiBarcode::new(0, vec![(0.0, 5.0), (2.0, 8.0)]);
        let curve = BettiCurve::from_barcode(&bc, 0.0, 10.0, 100);
        assert!(curve.value_at(3.0) >= 2);
        assert_eq!(curve.value_at(9.0), 0);
    }

    #[test]
    fn test_betti_summary() {
        let bc = BettiBarcode::new(0, vec![(0.0, 3.0), (1.0, 2.0)]);
        let summary = BettiSummary::from_barcode(&bc);
        assert!(summary.max >= 2);
        assert!(summary.total_persistence > 0.0);
        assert!(summary.mean > 0.0);
    }

    #[test]
    fn test_betti_summary_empty() {
        let bc = BettiBarcode::new(0, vec![]);
        let summary = BettiSummary::from_barcode(&bc);
        assert_eq!(summary.max, 0);
        assert_eq!(summary.mean, 0.0);
    }

    #[test]
    fn test_euler_curve() {
        let bc0 = BettiBarcode::new(0, vec![(0.0, 5.0), (0.0, 5.0)]);
        let bc1 = BettiBarcode::new(1, vec![(1.0, 3.0)]);
        let euler = EulerCurve::from_barcodes(&[bc0, bc1], 0.0, 5.0, 50);
        // At t=0: β₀=2, β₁=0 → χ=2
        assert!(euler.max() >= 2);
    }

    #[test]
    fn test_euler_curve_signs() {
        let bc0 = BettiBarcode::new(0, vec![(0.0, 10.0)]);
        let bc1 = BettiBarcode::new(1, vec![(0.0, 10.0)]);
        let bc2 = BettiBarcode::new(2, vec![(0.0, 10.0)]);
        let euler = EulerCurve::from_barcodes(&[bc0, bc1, bc2], 0.0, 10.0, 10);
        // χ = 1 - 1 + 1 = 1
        assert_eq!(euler.value_at_index(5), 1);
    }

    #[test]
    fn test_persistence_entropy_single_bar() {
        let bc = BettiBarcode::new(0, vec![(0.0, 1.0)]);
        let e = PersistenceEntropy::global_entropy(&bc);
        assert!(e.abs() < 1e-10); // single bar → entropy 0
    }

    #[test]
    fn test_persistence_entropy_equal_bars() {
        let bc = BettiBarcode::new(0, vec![(0.0, 1.0), (0.0, 1.0)]);
        let e = PersistenceEntropy::global_entropy(&bc);
        assert!((e - 2.0_f64.ln()).abs() < 1e-10);
    }

    #[test]
    fn test_persistence_entropy_empty() {
        let bc = BettiBarcode::new(0, vec![]);
        assert_eq!(PersistenceEntropy::global_entropy(&bc), 0.0);
    }
}
