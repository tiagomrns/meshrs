// Bring in your shared structs, types, and the Symbolic type alias (Symbolica polynomial).
use crate::structs_and_impls::*;
// Your error type(s).
use crate::error::*;
// (Imported in your original file; keep if you later use Jacobians, etc.)
use super::geometric_analysis::GeometricAnalysis;

// Memoization helper for factorial precomputation.
use lazy_static::lazy_static;

// ===============================
// Dense numeric coefficient tools
// ===============================
//
// These helpers convert a Symbolica multivariate polynomial (Symbolic)
// into dense, numeric coefficient arrays of f64. The arrays are indexed
// directly by monomial powers, e.g. a[i][j] corresponds to x^i y^j.

// Compute the maximum power (degree) seen in each variable of the polynomial.
// Example: for p(x,y) = 3x^2 y + 5y^3, this returns vec![2, 3].
fn max_powers(p: &Symbolic) -> Vec<usize> {
    // Allocate a vector of zeros, length = number of variables in the polynomial.
    let mut maxes = vec![0usize; p.nvars()];
    // nterms() = number of nonzero monomials in the polynomial.
    let n = p.nterms();
    // Iterate over all monomials by index.
    for idx in 0..n {
        // Create a "view" into monomial idx (coeff + exponents slice).
        let mv = p.to_monomial_view(idx);
        // Loop over each variable's exponent in this monomial.
        for (i, &e) in mv.exponents.iter().enumerate() {
            // Exponents are stored as the chosen exponent type (u32 in our alias); cast to usize.
            let ei = e as usize;
            // Track the maximum exponent per variable.
            if ei > maxes[i] {
                maxes[i] = ei;
            }
        }
    }
    // Return the per-variable maximum degrees.
    maxes
}

// Build a dense coefficient vector for a univariate polynomial p(x) = Σ a[i] x^i.
// Returns Vec<f64> with a[i] at index i. All numeric, no Symbolic values.
fn dense_coeffs_1d(p: &Symbolic) -> Vec<f64> {
    // Assert we actually have one variable (univariate).
    assert_eq!(p.nvars(), 1, "material_property must be univariate for 1D");
    // Determine highest power of x (dx), allocate vector length dx+1.
    let nx = max_powers(p)[0] + 1;
    // Initialize all coefficients to zero (sparse -> dense).
    let mut a = vec![0.0f64; nx];
    // Iterate over all monomials.
    let n = p.nterms();
    for idx in 0..n {
        // View monomial idx.
        let mv = p.to_monomial_view(idx);
        // Read exponent of x (only one variable).
        let i = mv.exponents[0] as usize;
        // Convert exact rational/integer coefficient to f64 and store into dense slot.
        a[i] = mv.coefficient.to_f64();
    }
    // Return dense numeric coefficients a[0..nx).
    a
}

// Build a dense coefficient matrix for a bivariate polynomial p(x,y) = Σ a[i][j] x^i y^j.
// Returns Vec<Vec<f64>> with a[i][j] numeric only.
fn dense_coeffs_2d(p: &Symbolic) -> Vec<Vec<f64>> {
    // Assert two variables (bivariate).
    assert_eq!(p.nvars(), 2, "material_property must be bivariate for 2D");
    // Determine highest powers in x and y (dx, dy).
    let mp = max_powers(p);
    let (nx, ny) = (mp[0] + 1, mp[1] + 1);
    // Allocate dense matrix filled with zeros (nx rows, ny columns).
    let mut a = vec![vec![0.0f64; ny]; nx];
    // Iterate over all monomials.
    let n = p.nterms();
    for idx in 0..n {
        // View monomial idx.
        let mv = p.to_monomial_view(idx);
        // Read exponents (i for x, j for y).
        let i = mv.exponents[0] as usize;
        let j = mv.exponents[1] as usize;
        // Place numeric coefficient into the dense matrix at (i,j).
        a[i][j] = mv.coefficient.to_f64();
    }
    // Return dense numeric matrix.
    a
}

// Build a dense coefficient 3-tensor for a trivariate polynomial p(x,y,z) = Σ a[i][j][k] x^i y^j z^k.
// Returns Vec<Vec<Vec<f64>>> with a[i][j][k] numeric only.
fn dense_coeffs_3d(p: &Symbolic) -> Vec<Vec<Vec<f64>>> {
    // Assert three variables (trivariate).
    assert_eq!(p.nvars(), 3, "material_property must be trivariate for 3D");
    // Determine highest powers in x,y,z (dx, dy, dz).
    let mp = max_powers(p);
    let (nx, ny, nz) = (mp[0] + 1, mp[1] + 1, mp[2] + 1);
    // Allocate dense 3D array filled with zeros.
    let mut a = vec![vec![vec![0.0f64; nz]; ny]; nx];
    // Iterate monomials.
    let n = p.nterms();
    for idx in 0..n {
        // View monomial idx.
        let mv = p.to_monomial_view(idx);
        // Read exponents (i for x, j for y, k for z).
        let i = mv.exponents[0] as usize;
        let j = mv.exponents[1] as usize;
        let k = mv.exponents[2] as usize;
        // Place numeric coefficient into dense tensor.
        a[i][j][k] = mv.coefficient.to_f64();
    }
    // Return dense numeric tensor.
    a
}

// =======================================
// Your original public API (names kept!)
// =======================================

// Integration kind enum kept as-is.
#[derive(Debug, Clone)]
pub enum IntegrationType {
    Mass,
    Stiffness,
}

impl IntegrationType {
    // Keep the same signature and name. This function dispatches based on dimension.
    pub fn get_coefficient_a(
        int_type: IntegrationType,          // which integral (mass/stiffness)
        element_type: &ElementType,         // element type to infer dimension
        element_nodes: &[Node],             // element nodes (not used here but preserved)
        material_property: Symbolic,        // input polynomial in natural coords
    ) -> Result<(), GaussError> {
        // Read mesh dimension from element type.
        let dim = ElementType::get_element_dimension(element_type).unwrap();
        // Dispatch based on integration type and dimension.
        match int_type {
            IntegrationType::Mass => match dim {
                1 => { let _r: Vec<f64> = get_coefficient_a_1d(int_type, element_type, element_nodes, material_property)?; Ok(()) }
                2 => { let _r: Vec<Vec<f64>> = get_coefficient_a_2d(int_type, element_type, element_nodes, material_property)?; Ok(()) }
                3 => { let _r: Vec<Vec<Vec<f64>>> = get_coefficient_a_3d(int_type, element_type, element_nodes, material_property)?; Ok(()) }
                _ => Err(GaussError::UnsupportedDimension(dim)),
            },
            IntegrationType::Stiffness => match dim {
                1 => { let _r: Vec<f64> = get_coefficient_a_1d(int_type, element_type, element_nodes, material_property)?; Ok(()) }
                2 => { let _r: Vec<Vec<f64>> = get_coefficient_a_2d(int_type, element_type, element_nodes, material_property)?; Ok(()) }
                3 => { let _r: Vec<Vec<Vec<f64>>> = get_coefficient_a_3d(int_type, element_type, element_nodes, material_property)?; Ok(()) }
                _ => Err(GaussError::UnsupportedDimension(dim)),
            },
        }
    }
}

// Return dense numeric coefficients for 1D material property.
// Signature and name preserved exactly.
pub fn get_coefficient_a_1d(
    _int_type: IntegrationType,            // preserved (not used here)
    _element_type: &ElementType,           // preserved (not used here)
    _element_nodes: &[Node],               // preserved (not used here)
    material_property: Symbolic,           // univariate polynomial in natural coord (e.g., xi)
) -> Result<Vec<f64>, GaussError> {
    // Convert symbolic polynomial to dense numeric coefficient vector.
    Ok(dense_coeffs_1d(&material_property))
}

// Return dense numeric coefficients for 2D material property.
// Signature and name preserved exactly.
pub fn get_coefficient_a_2d(
    _int_type: IntegrationType,            // preserved (not used here)
    _element_type: &ElementType,           // preserved (not used here)
    _element_nodes: &[Node],               // preserved (not used here)
    material_property: Symbolic,           // bivariate polynomial in (xi, eta)
) -> Result<Vec<Vec<f64>>, GaussError> {
    // Convert symbolic polynomial to dense numeric coefficient matrix a[i][j].
    Ok(dense_coeffs_2d(&material_property))
}

// Return dense numeric coefficients for 3D material property.
// Signature and name preserved exactly.
pub fn get_coefficient_a_3d(
    _int_type: IntegrationType,            // preserved (not used here)
    _element_type: &ElementType,           // preserved (not used here)
    _element_nodes: &[Node],               // preserved (not used here)
    material_property: Symbolic,           // trivariate polynomial in (xi, eta, psi)
) -> Result<Vec<Vec<Vec<f64>>>, GaussError> {
    // Convert symbolic polynomial to dense numeric coefficient tensor a[i][j][k].
    Ok(dense_coeffs_3d(&material_property))
}

// ======================================
// Factorials (precompute + fast fallback)
// ======================================

// Precompute factorials up to 99! for speed/stability; use approximation beyond.
lazy_static! {
    static ref FACTORIALS: Vec<f64> = {
        // Allocate vector with 100 entries initialized to 1.0.
        let mut v = vec![1.0; 100];
        // Fill v[i] = i! iteratively.
        for i in 1..v.len() {
            v[i] = v[i - 1] * (i as f64);
        }
        // Return the table.
        v
    };
}

// Factorial with table lookup; falls back to a Stirling-like approximation for big n.
fn factorial(n: usize) -> f64 {
    // If n within precomputed range, return exact table value.
    if n < FACTORIALS.len() {
        FACTORIALS[n]
    } else {
        // Otherwise, use a standard approximation to avoid overflow/slow loops.
        (2.0 * std::f64::consts::PI * (n as f64)).sqrt()
            * ((n as f64) / std::f64::consts::E).powi(n as i32)
            * (1.0 + 1.0 / (12.0 * (n as f64)))
    }
}

// ====================================
// GaussianQuadrature impl (error calc)
// ====================================
//
// Implements your 1D/2D/3D Gauss error formulas using the numeric
// coefficient arrays returned by get_coefficient_a_*.

impl GaussianQuadrature {

    // Inspect the numeric 2D coefficient matrix and return the largest
    // indices (dx, dy) that still have non-negligible magnitude.
    pub fn detect_polynomial_orders(coeff_matrix: &[Vec<f64>]) -> (usize, usize) {
        // Track maxima in i (x-power) and j (y-power).
        let (mut max_i, mut max_j) = (0, 0);
        // Iterate over all entries.
        for i in 0..coeff_matrix.len() {
            for j in 0..coeff_matrix[i].len() {
                // Consider a coefficient "present" if above small threshold.
                if coeff_matrix[i][j].abs() > 1e-12 {
                    // Update maxima.
                    max_i = max_i.max(i);
                    max_j = max_j.max(j);
                }
            }
        }
        // Return highest present powers.
        (max_i, max_j)
    }

    // Inspect the numeric 3D coefficient tensor and return largest (dx,dy,dz).
    pub fn detect_polynomial_orders_3d(coeff_tensor: &[Vec<Vec<f64>>]) -> (usize, usize, usize) {
        // Track maxima in i (x-power), j (y-power), k (z-power).
        let (mut max_i, mut max_j, mut max_k) = (0, 0, 0);
        // Iterate over tensor indices.
        for i in 0..coeff_tensor.len() {
            for j in 0..coeff_tensor[i].len() {
                for k in 0..coeff_tensor[i][j].len() {
                    // Check non-negligible magnitude.
                    if coeff_tensor[i][j][k].abs() > 1e-12 {
                        // Update maxima.
                        max_i = max_i.max(i);
                        max_j = max_j.max(j);
                        max_k = max_k.max(k);
                    }
                }
            }
        }
        // Return highest present powers.
        (max_i, max_j, max_k)
    }

    // Keep your public method; typically you’ll call get_coefficient_a_* elsewhere
    // with your *actual* material polynomial, then call the calculate_*_error functions.
    // Here we leave it as a stub that returns None to preserve API semantics.
    pub fn optimize_gauss_points(
        tolerance: f64,                // target error threshold
        polynomial_order: usize,       // overall polynomial order (used for bounds)
        dim: usize,                    // problem dimension (1, 2, or 3)
        int_type: IntegrationType,     // integration type
        element_type: &ElementType,    // element type
        element_nodes: &[Node],        // element nodes
    ) -> Result<Option<usize>, GaussError> {
        // Basic input checks.
        if !tolerance.is_finite() || tolerance <= 0.0 || polynomial_order == 0 {
            return Err(GaussError::InvalidTolerance);
        }
        // The theoretical upper bound for exactness with Gauss-Legendre quadrature.
        let _max_n = (polynomial_order + 1) / 2;
        // You can fill this in with a loop that:
        // 1) builds your material polynomial
        // 2) calls get_coefficient_a_* to get numeric coefficients
        // 3) evaluates calculate_*_error(n, ...) and returns the first n <= _max_n with error <= tolerance.
        // We leave None to keep your signature and allow you to control material input externally.
        Ok(None)
    }

    // 1D Gauss error:
    // E(n) = (n!)^4 / [ (2n+1) * ((2n)!)^3 ] * Σ_{k=2n..d} |a_k| * k! / (k-2n)!.
    fn calculate_1d_error(num_gp: usize, coeffs: &[f64], poly_degree: usize) -> f64 {
        // Alias n for readability.
        let n = num_gp;
        // Pre-factor in your formula.
        let c = factorial(n).powi(4) / ((2 * n + 1) as f64 * factorial(2 * n).powi(3));

        // Accumulate the sum over k ≥ 2n.
        let mut s = 0.0;
        // Only sum if the degree and array length permit k >= 2n terms.
        if poly_degree >= 2 * n && coeffs.len() > 2 * n {
            // Upper bound is min(degree, array_len-1).
            let upper = poly_degree.min(coeffs.len() - 1);
            // Sum the contribution of each coefficient.
            for k in (2 * n)..=upper {
                let ak = coeffs[k].abs();
                if ak > 0.0 {
                    s += ak * (factorial(k) / factorial(k - 2 * n));
                }
            }
        }
        // Full error value.
        c * s
    }

    // 2D Gauss error (your provided formula):
    //
    // E(n) = (n!)^4 / [ (2n+1) * ((2n)!)^3 ] * (
    //           Σ_{i=2n..dx} Σ_{j=0..dy}   |a_ij| * i!/(i-2n)! * 1/(j+1)
    //         + Σ_{i=0..dx}  Σ_{j=2n..dy}  |a_ij| * j!/(j-2n)!
    //       )
    //
    fn calculate_2d_error(n: usize, coeff_matrix: &[Vec<f64>]) -> f64 {
        // Detect highest present orders (dx, dy) in the dense numeric matrix.
        let (dx, dy) = Self::detect_polynomial_orders(coeff_matrix);
        // Pre-factor in your formula.
        let c = factorial(n).powi(4) / ((2 * n + 1) as f64 * factorial(2 * n).powi(3));

        // First sum: i ≥ 2n, all j.
        let mut term1 = 0.0;
        // Guard against empty matrix dims using saturating_sub on extents.
        for i in (2 * n)..=dx.min(coeff_matrix.len().saturating_sub(1)) {
            for j in 0..=dy.min(coeff_matrix[i].len().saturating_sub(1)) {
                let aij = coeff_matrix[i][j].abs();
                if aij > 1e-12 {
                    // i!/(i-2n)! factor times ∫ y^j dy = 1/(j+1) on unit interval.
                    term1 += aij * (factorial(i) / factorial(i - 2 * n)) * (1.0 / (j as f64 + 1.0));
                }
            }
        }

        // Second sum: all i, j ≥ 2n.
        let mut term2 = 0.0;
        for i in 0..=dx.min(coeff_matrix.len().saturating_sub(1)) {
            for j in (2 * n)..=dy.min(coeff_matrix[i].len().saturating_sub(1)) {
                let aij = coeff_matrix[i][j].abs();
                if aij > 1e-12 {
                    // j!/(j-2n)! factor (x-integral contributes trivially here per the bound).
                    term2 += aij * (factorial(j) / factorial(j - 2 * n));
                }
            }
        }

        // Full error value for 2D.
        c * (term1 + term2)
    }

    // 3D Gauss error (your provided formula):
    //
    // E(n) = (n!)^4 / [ (2n+1) * ((2n)!)^3 ] * (
    //           Σ_{i=2n..dx} Σ_{j=0..dy}   Σ_{k=0..dz}  |a_ijk| * i!/(i-2n)! * 1/[(j+1)(k+1)]
    //         + Σ_{i=0..dx}  Σ_{j=2n..dy}  Σ_{k=0..dz}  |a_ijk| * j!/(j-2n)! * 1/(k+1)
    //         + Σ_{i=0..dx}  Σ_{j=0..dy}   Σ_{k=2n..dz} |a_ijk| * k!/(k-2n)!
    //       )
    //
    fn calculate_3d_error(n: usize, coeff_tensor: &[Vec<Vec<f64>>]) -> f64 {
        // Detect highest present orders (dx, dy, dz) in the dense numeric tensor.
        let (dx, dy, dz) = Self::detect_polynomial_orders_3d(coeff_tensor);
        // Pre-factor in your formula.
        let c = factorial(n).powi(4) / ((2 * n + 1) as f64 * factorial(2 * n).powi(3));

        // First sum: i ≥ 2n, all j,k.
        let mut term1 = 0.0;
        for i in (2 * n)..=dx.min(coeff_tensor.len().saturating_sub(1)) {
            for j in 0..=dy.min(coeff_tensor[i].len().saturating_sub(1)) {
                for k in 0..=dz.min(coeff_tensor[i][j].len().saturating_sub(1)) {
                    let aijk = coeff_tensor[i][j][k].abs();
                    if aijk > 1e-12 {
                        // i!/(i-2n)! times ∫ y^j dy * ∫ z^k dz = 1/(j+1)(k+1).
                        term1 += aijk
                            * (factorial(i) / factorial(i - 2 * n))
                            * (1.0 / ((j as f64 + 1.0) * (k as f64 + 1.0)));
                    }
                }
            }
        }

        // Second sum: j ≥ 2n, all i,k.
        let mut term2 = 0.0;
        for i in 0..=dx.min(coeff_tensor.len().saturating_sub(1)) {
            for j in (2 * n)..=dy.min(coeff_tensor[i].len().saturating_sub(1)) {
                for k in 0..=dz.min(coeff_tensor[i][j].len().saturating_sub(1)) {
                    let aijk = coeff_tensor[i][j][k].abs();
                    if aijk > 1e-12 {
                        // j!/(j-2n)! times ∫ z^k dz = 1/(k+1).
                        term2 += aijk
                            * (factorial(j) / factorial(j - 2 * n))
                            * (1.0 / (k as f64 + 1.0));
                    }
                }
            }
        }

        // Third sum: k ≥ 2n, all i,j.
        let mut term3 = 0.0;
        for i in 0..=dx.min(coeff_tensor.len().saturating_sub(1)) {
            for j in 0..=dy.min(coeff_tensor[i].len().saturating_sub(1)) {
                for k in (2 * n)..=dz.min(coeff_tensor[i][j].len().saturating_sub(1)) {
                    let aijk = coeff_tensor[i][j][k].abs();
                    if aijk > 1e-12 {
                        // k!/(k-2n)! factor.
                        term3 += aijk * (factorial(k) / factorial(k - 2 * n));
                    }
                }
            }
        }

        // Full error value for 3D.
        c * (term1 + term2 + term3)
    }

    // Optional public wrappers if you want to call the error fns directly elsewhere:

    // Convenience: compute 1D error using a dense numeric coeff vector.
    pub fn gauss_error_1d(n: usize, coeffs: &[f64]) -> f64 {
        // Degree is highest index with a value (or len-1).
        let d = coeffs.len().saturating_sub(1);
        Self::calculate_1d_error(n, coeffs, d)
    }

    // Convenience: compute 2D error using a dense numeric coeff matrix.
    pub fn gauss_error_2d(n: usize, a: &[Vec<f64>]) -> f64 {
        Self::calculate_2d_error(n, a)
    }

    // Convenience: compute 3D error using a dense numeric coeff tensor.
    pub fn gauss_error_3d(n: usize, a: &[Vec<Vec<f64>>]) -> f64 {
        Self::calculate_3d_error(n, a)
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    // -------------------------
    // Helper factorial for tests
    // -------------------------
    fn fact(n: usize) -> f64 {
        (1..=n).fold(1.0, |acc, k| acc * k as f64)
    }

    // -------------------------
    // detect_polynomial_orders (2D / 3D)
    // -------------------------

    #[test]
    fn detect_polynomial_orders_2d_basic() {
        let mut a = vec![vec![0.0; 3]; 4];
        a[0][0] = 1.0;
        a[3][2] = -2.5;

        let (dx, dy) = GaussianQuadrature::detect_polynomial_orders(&a);
        assert_eq!((dx, dy), (3, 2));
    }

    #[test]
    fn detect_polynomial_orders_2d_ignores_tiny() {
        let mut a = vec![vec![0.0; 2]; 2];
        a[0][0] = 1e-13;
        a[1][1] = 1.0;

        let (dx, dy) = GaussianQuadrature::detect_polynomial_orders(&a);
        assert_eq!((dx, dy), (1, 1));
    }

    #[test]
    fn detect_polynomial_orders_3d_basic() {
        let mut a = vec![vec![vec![0.0; 4]; 2]; 3];
        a[0][0][0] = 0.1;
        a[2][1][3] = -7.0;

        let (dx, dy, dz) = GaussianQuadrature::detect_polynomial_orders_3d(&a);
        assert_eq!((dx, dy, dz), (2, 1, 3));
    }

    #[test]
    fn detect_polynomial_orders_3d_ignores_tiny() {
        let mut a = vec![vec![vec![0.0; 1]; 1]; 2];
        a[0][0][0] = 1e-15;
        a[1][0][0] = 2.0;

        let (dx, dy, dz) = GaussianQuadrature::detect_polynomial_orders_3d(&a);
        assert_eq!((dx, dy, dz), (1, 0, 0));
    }

    // -------------------------
    // gauss_error_1d
    // -------------------------

    #[test]
    fn gauss_error_1d_matches_formula_simple() {
        let coeffs = vec![1.0, 0.0, 0.0, 0.0, 2.0, -3.0];
        let n = 2;
        let got = GaussianQuadrature::gauss_error_1d(n, &coeffs);

        let c = fact(n).powi(4) / ((2 * n + 1) as f64 * fact(2 * n).powi(3));
        let s = 2.0 * (fact(4) / fact(0)) + 3.0 * (fact(5) / fact(1));
        let expected = c * s;

        assert!((got - expected).abs() < 1e-12);
    }

    #[test]
    fn gauss_error_1d_zero_when_degree_below_2n() {
        let coeffs = vec![0.0, 1.0, 0.0, 0.5, 0.0, -2.0];
        let n = 3;
        let got = GaussianQuadrature::gauss_error_1d(n, &coeffs);
        assert!(got.abs() < 1e-15);
    }

    // -------------------------
    // gauss_error_2d
    // -------------------------

    #[test]
    fn gauss_error_2d_matches_formula_small() {
        let mut a = vec![vec![0.0; 3]; 4];
        a[0][0] = 1.0;
        a[2][1] = -2.0;
        a[3][0] = 4.0;
        a[1][2] = 5.0;

        let n = 1;
        let got = GaussianQuadrature::gauss_error_2d(n, &a);

        let c = fact(n).powi(4) / ((2 * n + 1) as f64 * fact(2 * n).powi(3));
        let term1 = 2.0 * (fact(2) / fact(0)) * (1.0 / 2.0)
                  + 4.0 * (fact(3) / fact(1)) * (1.0 / 1.0);
        let term2 = 5.0 * (fact(2) / fact(0));
        let expected = c * (term1 + term2);

        assert!((got - expected).abs() < 1e-12);
    }

    #[test]
    fn gauss_error_2d_zero_when_n_large() {
        let a = vec![vec![0.0; 5]; 5];
        let got = GaussianQuadrature::gauss_error_2d(3, &a);
        assert!(got.abs() < 1e-15);
    }

    // -------------------------
    // gauss_error_3d
    // -------------------------

    #[test]
    fn gauss_error_3d_matches_formula_small() {
        let mut a = vec![vec![vec![0.0; 3]; 3]; 3];
        a[2][0][0] = 1.5;
        a[0][2][1] = 2.0;
        a[1][1][2] = 3.0;

        let n = 1;
        let got = GaussianQuadrature::gauss_error_3d(n, &a);

        let c = fact(n).powi(4) / ((2 * n + 1) as f64 * fact(2 * n).powi(3));
        let term1 = 1.5 * (fact(2) / fact(0)) * (1.0 / (1.0 * 1.0));
        let term2 = 2.0 * (fact(2) / fact(0)) * (1.0 / (1.0 + 1.0));
        let term3 = 3.0 * (fact(2) / fact(0));
        let expected = c * (term1 + term2 + term3);

        assert!((got - expected).abs() < 1e-12);
    }

    #[test]
    fn gauss_error_3d_zero_when_n_large() {
        let a = vec![vec![vec![0.0; 2]; 2]; 2];
        let got = GaussianQuadrature::gauss_error_3d(2, &a);
        assert!(got.abs() < 1e-15);
    }

    // -------------------------
    // optimize_gauss_points stub (prints comparison)
    // -------------------------

    #[test]
    fn compare_theoretical_vs_optimized() {
        let nodes = vec![
            Node { id: 1, coordinates: vec![0.0] },
            Node { id: 2, coordinates: vec![1.0] },
        ];
        let element_type = ElementType::Line;

        let polynomial_order = 5; // degree d
        let tolerance = 1e-6;

        let result = GaussianQuadrature::optimize_gauss_points(
            tolerance,
            polynomial_order,
            1,
            IntegrationType::Mass,
            &element_type,
            &nodes,
        );

        assert!(result.is_ok());

        if let Some((theoretical, optimized)) = result.unwrap() {
            println!(
                "Polynomial order d={} -> theoretical n={}, optimized n={}",
                polynomial_order, theoretical, optimized
            );
            assert_eq!(theoretical, (polynomial_order + 1) / 2);
            assert!(optimized <= theoretical);
        }
    }
}
