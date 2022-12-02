//! This module performs polynomial operations on a mathematical field 𝓕.
//!
//! 𝓕 has cardinality 256 so there is a bijection between a byte and a field element.
//!
//! A polynomial 𝑃, of degree 𝑛, is fully defined by 𝑛+1 coefficients 𝑎ᵢ and
//! is written as: 𝑎ₙ𝑋ⁿ + 𝑎ₙ₋₁𝑋ⁿ⁻¹ + ... + 𝑎₂𝑋² + 𝑎₁𝑋 + 𝑎₀
//!
//! By convention, the degree of polynomial 𝑃 = 0 is −∞.

// ----------------
extern crate test;
use crate::galois_field_256::Element as Elt;

/// A polynomial 𝑃 = 𝑎ₙ𝑋ⁿ + 𝑎ₙ₋₁𝑋ⁿ⁻¹ + ... + 𝑎₂𝑋² + 𝑎₁𝑋 + 𝑎₀ is represented as the vector
/// `[𝑎₀, ..., 𝑎ₙ]`.
///
/// A polynomial of degree 𝑛 is therefore represented by 𝑛+1 coefficients.
///
/// The last coefficient in the vector (𝑎ₙ) is never equal to zero, otherwise the degree of the
/// polynomial would be less than 𝑛.
///
/// The zero polynomial is represented by the empty vector.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct GF256Polynomial(Vec<Elt>);

impl std::ops::AddAssign<&Self> for GF256Polynomial {
    /// Performs in-place addition with another, borrowed, polynomial
    ///
    /// Addition uses the addition of 𝓕 on coefficients of the same degree
    ///
    /// # Parameters
    ///
    /// | Parameter | Description | Notes |
    /// | --------- | ----------- | ----- |
    /// | `rhs` | operand for addition | |
    fn add_assign(&mut self, rhs: &Self) {
        match (self.degree(), rhs.degree()) {
            (None, _) => {
                self.0 = rhs.0.clone();
            }
            (_, None) => {}
            (Some(self_degree), Some(rhs_degree)) => {
                if self_degree < rhs_degree {
                    self.0.resize(rhs_degree + 1, Elt::new(0));
                }

                for (idx, &coef) in rhs.coefficients().iter().enumerate() {
                    self.0[idx] = self.0[idx] + coef;
                }

                while !self.0.is_empty() && self.0[self.0.len() - 1].value() == 0 {
                    self.0.pop();
                }
            }
        }
    }
}

impl std::ops::AddAssign<Self> for GF256Polynomial {
    /// Performs in-place addition with another, owned, polynomial
    ///
    /// Addition uses the addition of 𝓕 on coefficients of the same degree
    ///
    /// Taking ownership of the operand makes it possible to reuse the allocated storage for the
    /// largest vector
    ///
    /// # Parameters
    ///
    /// | Parameter | Description | Notes |
    /// | --------- | ----------- | ----- |
    /// | `rhs` | operand for addition | |
    fn add_assign(&mut self, mut rhs: Self) {
        match (self.degree(), rhs.degree()) {
            (None, _) => {
                self.0 = rhs.0;
            }
            (_, None) => {}
            (Some(self_degree), Some(rhs_degree)) => {
                if self_degree < rhs_degree {
                    // owning the vector with largest capacity
                    std::mem::swap(&mut self.0, &mut rhs.0);
                }

                for (idx, &coef) in rhs.coefficients().iter().enumerate() {
                    self.0[idx] = self.0[idx] + coef;
                }

                while !self.0.is_empty() && self.0[self.0.len() - 1].value() == 0 {
                    self.0.pop();
                }
            }
        }
    }
}

impl<'a, 'b> std::ops::Add<&'b GF256Polynomial> for &'a GF256Polynomial {
    type Output = GF256Polynomial;
    /// Performs addition with another polynomial
    ///
    /// Addition uses the addition of 𝓕 on coefficients of the same degree
    ///
    /// # Parameters
    ///
    /// | Parameter | Description | Notes |
    /// | --------- | ----------- | ----- |
    /// | `rhs` | operand for addition | |
    fn add(self, rhs: &'b GF256Polynomial) -> GF256Polynomial {
        let mut ret: GF256Polynomial = self.clone();
        ret += rhs;
        ret
    }
}

impl<'a> std::ops::Neg for &'a GF256Polynomial {
    type Output = Self;
    /// Computes the opposite of a polynomial
    ///
    /// # Output
    ///
    /// A new polynomial that, added with `self`, would result in the null polynomial.
    fn neg(self) -> Self {
        // assuming GF256 elements are self-opposites
        self
    }
}

impl std::ops::SubAssign<&Self> for GF256Polynomial {
    /// Performs in-place subtraction with another, borrowed, polynomial
    ///
    /// Subtraction uses the subtraction of 𝓕 on coefficients of the same degree
    ///
    /// # Parameters
    ///
    /// | Parameter | Description | Notes |
    /// | --------- | ----------- | ----- |
    /// | `rhs` | operand for subtraction | |
    fn sub_assign(&mut self, rhs: &Self) {
        // Assuming GF256 addition is subtraction
        *self += -rhs;
    }
}

impl std::ops::SubAssign<Self> for GF256Polynomial {
    /// Performs in-place subtraction with another, owned, polynomial
    ///
    /// Taking ownership of the operand makes it possible to reuse the allocated storage for the
    /// largest vector
    ///
    /// Subtraction uses the subtraction of 𝓕 on coefficients of the same degree
    ///
    /// # Parameters
    ///
    /// | Parameter | Description | Notes |
    /// | --------- | ----------- | ----- |
    /// | `rhs` | operand for subtraction | |
    fn sub_assign(&mut self, rhs: Self) {
        // Assuming GF256 addition is subtraction
        *self += -&rhs;
    }
}

impl<'a, 'b> std::ops::Sub<&'b GF256Polynomial> for &'a GF256Polynomial {
    type Output = GF256Polynomial;
    /// Performs subtraction with another polynomial
    ///
    /// Subtraction uses the subtraction of 𝓕 on coefficients of the same degree
    ///
    /// # Parameters
    ///
    /// | Parameter | Description | Notes |
    /// | --------- | ----------- | ----- |
    /// | `rhs` | operand for subtraction | |
    fn sub(self, rhs: &'b GF256Polynomial) -> GF256Polynomial {
        let mut ret: GF256Polynomial = self.clone();
        ret -= rhs;
        ret
    }
}

impl std::ops::MulAssign<Elt> for GF256Polynomial {
    /// Multiplies in place each coefficient of the polynomial by a constant factor of 𝓕.
    ///
    /// # Parameters
    ///
    /// | Parameter | Description | Notes |
    /// | --------- | ----------- | ----- |
    /// | `factor` | factor for multiplication | |
    fn mul_assign(&mut self, factor: Elt) {
        if factor == Elt::ZERO {
            self.0.clear();
        } else {
            for coef in &mut self.0 {
                *coef = *coef * factor;
            }
        }
    }
}

impl std::ops::Mul<Elt> for &GF256Polynomial {
    type Output = GF256Polynomial;

    /// Multiplies each coefficiont of the polynomial by a constant factor of 𝓕.
    ///
    /// # Parameters
    ///
    /// | Parameter | Description | Notes |
    /// | --------- | ----------- | ----- |
    /// | `factor` | factor for multiplication | |
    ///
    /// # Output
    ///
    /// New polynomial, product of `factor` × `self`.
    fn mul(self, factor: Elt) -> GF256Polynomial {
        let mut result = self.clone();
        result *= factor;
        result
    }
}

impl<'a, 'b> std::ops::Mul<&'b GF256Polynomial> for &'a GF256Polynomial {
    type Output = GF256Polynomial;
    /// Performs multiplication with another polynomial.
    ///
    /// Multiplication of polynomials uses the addition and multiplication of 𝓕.
    ///
    /// # Parameters
    ///
    /// | Parameter | Description | Notes |
    /// | --------- | ----------- | ----- |
    /// | `rhs` | operand for multiplication | |
    ///
    /// # Output
    ///
    /// New polynomial, product of `self` and `rhs`.
    fn mul(self, rhs: &'b GF256Polynomial) -> GF256Polynomial {
        match (self.degree(), rhs.degree()) {
            (None, _) | (_, None) => GF256Polynomial::new(Vec::new()),
            (Some(self_degree), Some(rhs_degree)) => {
                let product_degree = self_degree + rhs_degree;
                let mut result = Vec::with_capacity(product_degree + 1);
                result.resize(product_degree + 1, Elt::new(0));

                for (degree, item) in result.iter_mut().enumerate().take(product_degree + 1) {
                    let mut coefficient = Elt::new(0);

                    for index in 0..=degree {
                        if index >= self.coefficients().len() {
                            continue;
                        }
                        if degree - index >= rhs.coefficients().len() {
                            continue;
                        }

                        coefficient = coefficient
                            + (self.coefficients()[index] * rhs.coefficients()[degree - index]);
                    }

                    *item = coefficient;
                }

                GF256Polynomial::new(result)
            }
        }
    }
}

impl std::ops::MulAssign<&Self> for GF256Polynomial {
    /// Multiplies self with another polynomial.
    ///
    /// Multiplication of polynomials uses the addition and multiplication of 𝓕.
    ///
    /// # Parameters
    ///
    /// | Parameter | Description | Notes |
    /// | --------- | ----------- | ----- |
    /// | `rhs` | operand for multiplication | |
    fn mul_assign(&mut self, rhs: &Self) {
        *self = &self.clone() * rhs;
    }
}

impl GF256Polynomial {
    /// Retrieves the vector of coefficients of the polynomial
    ///
    /// # Output
    ///
    /// * Vector `[𝑎₀, ..., 𝑎ₙ]` where 𝑛 is the degree of the polynomial.
    /// * Empty vector for the null polynomial.
    #[inline]
    const fn coefficients(&self) -> &Vec<Elt> {
        &self.0
    }

    /// Constructs a new polynomial from its coefficients
    ///
    /// # Parameters
    ///
    /// | Parameter | Description | Notes |
    /// | --------- | ----------- | ----- |
    /// | `coefficients` | Vector `[𝑎₀, ..., 𝑎ₙ]` where 𝑛 is the degree of the polynomial | can be empty for null polynomial |
    ///
    /// # Output
    ///
    /// Newly created polynomial
    #[inline]
    fn new(coefficients: Vec<Elt>) -> Self {
        Self(coefficients)
    }

    /// Creates a random polynomial of degree at most 𝑛, and, optionnally, known value when 𝑥 = 0.
    ///
    /// If you want to construct the null polynomial, use a maximal degree of zero and constant
    /// coefficient of 0.
    ///
    /// # Parameters
    ///
    /// | Parameter | Description | Notes |
    /// | --------- | ----------- | ----- |
    /// | `max_degree` | Value of 𝑛 (maximal allowable degree of the polynomial) | use 0 for the null polynomial |
    /// | `constant_coefficient` | Optionnally, desired byte-mapped value of the coefficient 𝑎₀ | use `Some(0)` for the null polynomial |
    /// | `rng` | A random number generator | |
    ///
    /// # Output
    ///
    /// Newly created polynomial
    pub fn new_random<RNG>(max_degree: u8, constant_coefficient: Option<u8>, rng: &mut RNG) -> Self
    where
        RNG: rand::Rng,
    {
        use rand::distributions::Distribution;

        match (max_degree, constant_coefficient) {
            (0, Some(elt)) if elt == 0 => Self::new(Vec::new()),
            (0, Some(val)) => Self::new(vec![Elt::new(val)]),
            (_, Some(val)) => Self::new(
                std::iter::once(Elt::new(val))
                    .chain(
                        rand::distributions::Uniform::new_inclusive(0_u8, 255_u8)
                            .sample_iter(&mut *rng)
                            .take(max_degree as usize)
                            .map(Elt::new),
                    )
                    .collect(),
            ),
            (_, None) => {
                let coefficients: Vec<Elt> =
                    rand::distributions::Uniform::new_inclusive(0_u8, 255_u8)
                        .sample_iter(&mut *rng)
                        .take(max_degree as usize + 1)
                        .map(Elt::new)
                        .collect();
                Self::new(coefficients)
            }
        }
    }

    /// Retrieves the degree 𝑛 of the polynomial
    ///
    /// # Output
    ///
    /// - `Some(𝑛)`: the polynomial is of degree 𝑛
    /// - `None`: the polynomial is null, of degree −∞.
    pub fn degree(&self) -> Option<usize> {
        if self.coefficients().is_empty() {
            None
        } else {
            Some(self.coefficients().len() - 1)
        }
    }

    /// Compute a polynomial of degree 𝑛 using at least 𝑛+1 points. Points are defined by a pair of
    /// bytes that each maps to a value of 𝓕.
    ///
    /// This method uses Lagrange interpolation. Given a set of points (𝑥ᵢ, 𝑦ᵢ), then:
    ///
    /// * The partial Lagrange polynoms 𝑃ᵢ are, for each 𝑖, the products of (𝑋 - 𝑥ⱼ)/(𝑥ᵢ - 𝑥ⱼ) for 𝑗 ≠ 𝑖.  
    ///   These polynoms have value 1 when 𝑥 = 𝑥ᵢ, and 0 when 𝑥 = 𝑥ⱼ, 𝑗 ≠ 𝑖
    /// * The reconstituted polynom is the sum of 𝑦ᵢ𝑃ᵢ
    ///
    /// # Parameters
    ///
    /// | Parameter | Description | Notes |
    /// | --------- | ----------- | ----- |
    /// | `x_values` | vector of at least 𝑛+1 abscissae (𝑥ᵢ) used for interpolation | Values 𝑥ᵢ must be distinct |
    /// | `y_values` | vector of ordinates, of same length as `x_values` | |
    /// | `degree` | value of 𝑛 | |
    ///
    /// # Output
    ///
    /// New polynomial satisfying 𝑃(𝑥ᵢ)= 𝑦ᵢ for each value of 𝑖
    pub fn new_from_points(x_values: &[u8], y_values: &[u8], degree: u8) -> Self {
        assert!(
            x_values.len() == y_values.len(),
            "x_values and y_values must have same length"
        );
        assert!(x_values.len() > degree.into(), "not enough points provided");

        if x_values.is_empty() {
            return Self::new(Vec::new());
        }

        let mut result = Self::new(Vec::with_capacity(degree as usize + 1));

        for i in 0_usize..=(degree as usize) {
            let mut partial_lagrange_polynomial =
                Self::new(Vec::with_capacity(degree as usize + 1));
            partial_lagrange_polynomial.0.push(Elt::ONE);

            // Compute 𝑃ᵢ
            for j in 0_usize..=(degree as usize) {
                if i == j {
                    continue;
                }

                assert!(
                    x_values[i] != x_values[j],
                    "Multiple points with same abscissa"
                );

                // 𝑋 - 𝑥ⱼ is represented as (−𝑥ⱼ, 1)
                let mut numerator = Self::new(vec![-Elt::new(x_values[j]), Elt::ONE]);
                numerator *= (Elt::new(x_values[i]) - Elt::new(x_values[j]))
                    .inverse()
                    .unwrap();

                // multiply by 1/(𝑥ᵢ - 𝑥ⱼ)
                partial_lagrange_polynomial *= &numerator;
            }

            partial_lagrange_polynomial *= Elt::new(y_values[i]);

            // Add to result, times 𝑦ᵢ
            result += &partial_lagrange_polynomial;
        }

        result
    }

    /// Computes the byte-mapping of 𝑃(𝑥) given the byte-mapping of 𝑥.
    ///
    /// # Parameters
    ///
    /// | Parameter | Description | Notes |
    /// | --------- | ----------- | ----- |
    /// | `x` | Value of 𝑥 | |
    ///
    /// # Output
    ///
    /// Value of 𝑃(𝑥)
    pub fn value_at(&self, x: u8) -> u8 {
        // Optimization for x=0 : take constant coefficient
        if x == Elt::ZERO.value() {
            match self.degree() {
                None => return 0,
                Some(_) => return self.coefficients()[0].value(),
            };
        }

        let mut result = Elt::new(0);
        for i in 0..self.coefficients().len() {
            result = result + (self.coefficients()[i] * Elt::new(x).power(i as u8));
        }

        result.value()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn polynomial_arithmetic_is_working() {
        let mut rng = rand::thread_rng();

        let p1 = GF256Polynomial::new_random(5, Some(10), &mut rng);
        let p2 = GF256Polynomial::new_random(3, Some(4), &mut rng);

        assert_eq!(*(&p1 + &p1).coefficients(), Vec::new());
        assert_eq!(*(&p1 - &p1).coefficients(), Vec::new());
        assert_eq!(p1.degree(), Some(5), "p1 {p1:?} degree is not 5");
        assert_ne!(p1.coefficients()[5], Elt::new(0));
        assert_eq!(p1.coefficients()[0], Elt::new(10));

        let p1_times_p2 = &p1 * &p2;

        assert_eq!(
            p1_times_p2.degree().unwrap(),
            p1.degree().unwrap() + p2.degree().unwrap(),
            "Degrees do not add up: p1 = {:?}, p2 = {:?}, product = {:?}",
            p1.coefficients(),
            p2.coefficients(),
            p1_times_p2.coefficients()
        );
        assert_eq!(
            p1_times_p2.coefficients()[0],
            p1.coefficients()[0] * p2.coefficients()[0]
        );
        assert_eq!(
            p1_times_p2.coefficients()[8],
            p1.coefficients()[5] * p2.coefficients()[3]
        );
        assert_eq!(
            p1_times_p2.coefficients()[1],
            (p1.coefficients()[0] * p2.coefficients()[1])
                + (p1.coefficients()[1] * p2.coefficients()[0])
        );
    }

    #[test]
    fn polynomial_interpolation_works_as_expected() {
        let mut rng = rand::thread_rng();

        for degree in 1..10 {
            let random_polynomial = GF256Polynomial::new_random(degree, None, &mut rng);
            let x_values: Vec<u8> = (1..=(degree + 1)).map(|x| x * 2).collect();
            let y_values: Vec<u8> = (1..=(degree + 1))
                .map(|x| random_polynomial.value_at(x * 2))
                .collect();

            let reconstituted_p =
                GF256Polynomial::new_from_points(&x_values, &y_values, x_values.len() as u8 - 1);
            assert_eq!(
                random_polynomial,
                reconstituted_p,
                "Cannot reconstitute polynom. original = {:?}, points = {:?}, reconstituted = {:?}",
                random_polynomial.coefficients(),
                x_values
                    .iter()
                    .zip(y_values.iter())
                    .flat_map(|(x, y)| [*x, *y])
                    .collect::<Vec<u8>>(),
                reconstituted_p.coefficients()
            );

            for (x, y) in x_values.iter().zip(y_values.iter()) {
                let recovered_value = reconstituted_p.value_at(*x);
                assert_eq!(recovered_value, *y);
            }
        }
    }

    #[bench]
    fn bench_polynomial_creation(b: &mut test::Bencher) {
        let mut rng = rand::thread_rng();
        b.iter(move || {
            for degree in 2..=255 {
                test::black_box(GF256Polynomial::new_random(degree, Some(0), &mut rng));
            }
        });
    }

    #[bench]
    fn bench_polynomial_interpolation(b: &mut test::Bencher) {
        let mut rng = rand::thread_rng();
        let degree = 30;
        let poly = GF256Polynomial::new_random(degree, Some(0), &mut rng);
        let x_values: Vec<u8> = (1..=(degree + 1)).collect();
        let y_values: Vec<u8> = x_values.iter().map(|x| poly.value_at(*x)).collect();
        b.iter(move || {
            test::black_box(GF256Polynomial::new_from_points(
                &x_values, &y_values, degree,
            ));
        });
    }
}
