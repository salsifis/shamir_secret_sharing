//! This module performs arithmetics in a mathematical field with 256 elements and provides a
//! mapping between a single byte and an element of that field.
//!
//! # Mathematical definition of a field
//!
//! A field 𝓕 is a set on which addition (+) and multiplication (×) are defined and satisfy the
//! following requirements:
//!
//! - 𝓕 is a ring under addition and multiplication:
//!   - 𝓕 is an abelian group (or commutative group) under addition:
//!     - 𝓕 is a group under addition:
//!       - 𝓕 is a monoid under addition:
//!         - 𝓕 is a semi-group under addition:
//!           - 𝓕 is a magma under addition, which means addition is defined for each pair of values of 𝓕
//!           - addition in 𝓕 is associative: ∀ (𝑥,𝑦,𝑧) ∈ 𝓕³, (𝑥 + 𝑦) + 𝑧 = 𝑥 + (𝑦 + 𝑧)
//!         - There is an additive identity element in 𝓕, noted 0, such that ∀ 𝑦 ∈ 𝓕, 0 + 𝑦 = 𝑦 + 0 = 𝑦
//!       - for every element 𝑥 in 𝓕, there is an additive inverse element 𝑦 in 𝓕 such that 𝑥 + 𝑦 = 𝑦 + 𝑥 = 0
//!     - addition in 𝓕 is commutative: ∀ (𝑥,𝑦) ∈ 𝓕², 𝑥 + 𝑦 = 𝑦 + 𝑥
//!   - 𝓕 is a monoid under multiplication:
//!     - 𝓕 is a semi-group under multiplication:
//!       - 𝓕 is a magma under multiplication, which means multiplication is defined for each pair of values of 𝓕
//!       - multiplication in 𝓕 is associative: ∀ (𝑥,𝑦,𝑧) ∈ 𝓕³, (𝑥 × 𝑦) × 𝑧 = 𝑥 × (𝑦 × 𝑧)
//!     - There is an multiplicative identity element in 𝓕, noted 1, such that ∀ 𝑦 ∈ 𝓕, 1 × 𝑦 = 𝑦 × 1 = 𝑦
//!   - multiplication is left-distributive and right-distributive with respect to addition:  
//!     ∀ (𝑥,𝑦,𝑧) ∈ 𝓕³, 𝑥 × (𝑦 + 𝑧) = (𝑥 × 𝑦) + (𝑥 × 𝑧) and (𝑥 + 𝑦) × 𝑧 = (𝑥 × 𝑧) + (𝑦 × 𝑧).
//! - 𝓕 ∖ {0} is an abelian group under multiplication, which means that, additionnally:
//!   - for every element 𝑥 in 𝓕 except 0, there is a multiplicative inverse element 𝑦 in 𝓕 such that 𝑥 × 𝑦 = 𝑦 × 𝑥 = 1
//!   - multiplication in 𝓕 is commutative
//!
//! # The Galois field GF(256)
//!
//! Each element in the Galois field can be thought of as a polynomial of degree up to 7.
//!
//! Coefficients of the polynomial pertain to group (ℤ/2ℤ, +), so they have two distinct possible
//! values, and there are exactly 256 possible polynomials as there are 8 coefficients total.
//!
//! ## Addition
//!
//! The addition table in ℤ/2ℤ is as follows:
//!
//! | 𝑥 + 𝑦 | 𝑦 : 0 | 𝑦 : 1 |
//! | ----- | ----- | ----- |
//! | 𝑥 : 0 | 0     | 1     |
//! | 𝑥 : 1 | 1     | 0     |
//!
//! As for any polynomial structures, addition between polynomials just adds coefficients of the
//! same degree.  
//! For example (𝑋⁷ + 𝑋⁵ + 𝑋 + 1) + (𝑋³ + 𝑋² + 𝑋) = (𝑋⁷ + 𝑋⁵ + 𝑋³ + 𝑋² + 1)
//!
//! Addition is commutative and associative. The polynomial 0 is the neutral element.
//!
//! ## Multiplication
//!
//! Multiplication uses the standard multiplication of polynomials:  
//! if 𝑎ᵢ are the coefficients of polynomial 𝑃 and 𝑏ᵢ the coefficients of polynomial 𝑄, then
//! the coefficient of degree 𝑛 of 𝑃𝑄 is the sum of 𝑎ᵢ×𝑏ⱼ such that 𝑖+𝑗=𝑛.
//!
//! The multiplication table in ℤ/2ℤ is as follows:
//!
//! | 𝑥 × 𝑦 | 𝑦 : 0 | 𝑦 : 1 |
//! | ----- | ----- | ----- |
//! | 𝑥 : 0 | 0     | 0     |
//! | 𝑥 : 1 | 0     | 1     |
//!
//! Since the degree of the product polynomial could exceed 7, we sutract multiples of
//! an irreductible polynomial of degree 8. A list of irreductible polynomials can be computed
//! using a function in this module. Here we use 𝑋⁸ + 𝑋⁷ + 𝑋⁶ + 𝑋³ + 𝑋² + 𝑋  + 1.
//!
//! Multiplication is commutative and associative. 1 is the neutral element.
//!
//! Tests ensure that for each non-zero element, there is exactly one inverse element, such that
//! the product of both elements is 1.
//!
//! For division of elements, the table of inverses is built at compile time, and division is
//! redirected to a multiplication with the inverse value.
//!
//! # About this implementation
//!
//! This implementation uses lookup tables that are evaluated at compilation time for operations
//! such as multiplication, powers and division. This is feasible for GF(256) as the lookup tables
//! are small.

// ----------------
extern crate test;

/// An element is conceptually a polynomial of degree ⩽ 7, with coefficients in ℤ/2ℤ (0 or 1).
///
/// Any of the 256 polynomials are mapped to a specific byte value, with each bit position
/// representing a coefficient of a specific degree.
///
/// For example, 163 (0b10100011) is used to represent the element that can be
/// thought of as the polynomial 𝑋⁷ + 𝑋⁵ + 𝑋 + 1.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct Element(u8);

impl std::ops::Add for Element {
    type Output = Self;

    /// Performs addition with another element.
    ///
    /// # Parameters
    ///
    /// | Parameter | Description | Notes |
    /// | --------- | ----------- | ----- |
    /// | `rhs` | operand for addition | |
    ///
    /// # Output
    ///
    /// New element, sum of `self` and `rhs`.
    #[allow(clippy::suspicious_arithmetic_impl)]
    fn add(self, rhs: Self) -> Self {
        // Addition of coefficients in ℤ/2ℤ is equivalent to XORing the values.
        Self::new(self.0 ^ rhs.0)
    }
}

impl std::ops::Neg for Element {
    type Output = Self;

    /// Computes the opposive (additive inverse) of an element
    ///
    /// # Output
    ///
    /// A new element that, added with `self`, would result in the neutral addition element
    fn neg(self) -> Self {
        // Elements in ℤ/2ℤ are their own opposite
        self
    }
}

impl std::ops::Sub for Element {
    type Output = Self;

    /// Performs subtraction with another element.
    ///
    /// # Parameters
    ///
    /// | Parameter | Description | Notes |
    /// | --------- | ----------- | ----- |
    /// | `rhs` | operand for subtraction | |
    ///
    /// # Output
    ///
    /// New element, difference of `self` and `rhs`.
    fn sub(self, rhs: Self) -> Self {
        self + -rhs
    }
}

impl std::ops::Mul for Element {
    type Output = Self;

    /// Computes the product with an element by direct lookup into a table
    ///
    /// # Parameters
    ///
    /// | Parameter | Description | Notes |
    /// | --------- | ----------- | ----- |
    /// | `rhs` | operand for multiplication | |
    ///
    /// # Output
    ///
    /// New element, product of `self` and `rhs`.
    fn mul(self, rhs: Self) -> Self {
        Self::ALL_PRODUCTS[self.value() as usize][rhs.value() as usize]
    }
}

impl std::ops::Div for Element {
    type Output = Option<Self>;

    /// Performs division by another element
    ///
    /// # Parameters
    ///
    /// | Parameter | Description | Notes |
    /// | --------- | ----------- | ----- |
    /// | `rhs` | operand for division | |
    ///
    /// # Output
    ///
    /// - If `rhs` is zero, then None
    /// - else, some new element, quotient of `self` by `rhs`.
    fn div(self, rhs: Self) -> Option<Self> {
        Self::ALL_ELEMENT_INVERSES[rhs.value() as usize]
            .map(|x| Self::ALL_PRODUCTS[self.value() as usize][x.value() as usize])
    }
}

impl Element {
    /// Representation of a polynomial of degree 8, which is irreductible in the set of
    /// polynomials with coefficients among 0 and 1.
    ///
    /// This means that no other polynomial of lesser degree divides it except for polynomials of degree 0.
    ///
    /// For example, an even number like 0x1c4, representing the polynomial 𝑋⁸ + 𝑋⁷ + 𝑋⁶ + 𝑋², is
    /// divisible by 0x1 (𝑋), which is a polynomial of degree 1, without remainder.
    ///
    /// A polynomial like 0x111 (𝑋⁸ + 𝑋⁴ + 1) is equal to (𝑋² + 𝑋 + 1)(𝑋⁶ + 𝑋⁵ + 𝑋³ + 𝑋 + 1) and hence is
    /// not irreductible (given the addition operation here).
    ///
    /// the function `compute_irreductible_polynomials()` lists these polynomials as satisfying this
    /// requirement:
    ///
    /// | Polynomial                      | Binary          | Number |
    /// | ------------------------------- | --------------- | ------ |
    /// | 𝑋⁸ + 𝑋⁴ + 𝑋³ + 𝑋  + 1           | `0b1_0001_1011` | 0x11b  |
    /// | 𝑋⁸ + 𝑋⁴ + 𝑋³ + 𝑋² + 1           | `0b1_0001_1101` | 0x11d  |
    /// | 𝑋⁸ + 𝑋⁵ + 𝑋³ + 𝑋  + 1           | `0b1_0010_1011` | 0x12b  |
    /// | 𝑋⁸ + 𝑋⁵ + 𝑋³ + 𝑋² + 1           | `0b1_0010_1101` | 0x12d  |
    /// | 𝑋⁸ + 𝑋⁵ + 𝑋⁴ + 𝑋³ + 1           | `0b1_0011_1001` | 0x139  |
    /// | 𝑋⁸ + 𝑋⁵ + 𝑋⁴ + 𝑋³ + 𝑋² + 𝑋  + 1 | `0b1_0011_1111` | 0x13f  |
    /// | 𝑋⁸ + 𝑋⁶ + 𝑋³ + 𝑋² + 1           | `0b1_0100_1101` | 0x14d  |
    /// | 𝑋⁸ + 𝑋⁶ + 𝑋⁴ + 𝑋³ + 𝑋² + 𝑋  + 1 | `0b1_0101_1111` | 0x15f  |
    /// | 𝑋⁸ + 𝑋⁶ + 𝑋⁵ + 𝑋  + 1           | `0b1_0110_0011` | 0x163  |
    /// | 𝑋⁸ + 𝑋⁶ + 𝑋⁵ + 𝑋² + 1           | `0b1_0110_0101` | 0x165  |
    /// | 𝑋⁸ + 𝑋⁶ + 𝑋⁵ + 𝑋³ + 1           | `0b1_0110_1001` | 0x169  |
    /// | 𝑋⁸ + 𝑋⁶ + 𝑋⁵ + 𝑋⁴ + 1           | `0b1_0111_0001` | 0x171  |
    /// | 𝑋⁸ + 𝑋⁶ + 𝑋⁵ + 𝑋⁴ + 𝑋² + 𝑋  + 1 | `0b1_0111_0111` | 0x177  |
    /// | 𝑋⁸ + 𝑋⁶ + 𝑋⁵ + 𝑋⁴ + 𝑋³ + 𝑋  + 1 | `0b1_0111_1011` | 0x17b  |
    /// | 𝑋⁸ + 𝑋⁷ + 𝑋² + 𝑋  + 1           | `0b1_1000_0111` | 0x187  |
    /// | 𝑋⁸ + 𝑋⁷ + 𝑋³ + 𝑋  + 1           | `0b1_1000_1011` | 0x18b  |
    /// | 𝑋⁸ + 𝑋⁷ + 𝑋³ + 𝑋² + 1           | `0b1_1000_1101` | 0x18d  |
    /// | 𝑋⁸ + 𝑋⁷ + 𝑋⁴ + 𝑋³ + 𝑋² + 𝑋  + 1 | `0b1_1001_1111` | 0x19f  |
    /// | 𝑋⁸ + 𝑋⁷ + 𝑋⁵ + 𝑋  + 1           | `0b1_1010_0011` | 0x1a3  |
    /// | 𝑋⁸ + 𝑋⁷ + 𝑋⁵ + 𝑋³ + 1           | `0b1_1010_1001` | 0x1a9  |
    /// | 𝑋⁸ + 𝑋⁷ + 𝑋⁵ + 𝑋⁴ + 1           | `0b1_1011_0001` | 0x1b1  |
    /// | 𝑋⁸ + 𝑋⁷ + 𝑋⁵ + 𝑋⁴ + 𝑋³ + 𝑋² + 1 | `0b1_1011_1101` | 0x1bd  |
    /// | 𝑋⁸ + 𝑋⁷ + 𝑋⁶ + 𝑋  + 1           | `0b1_1100_0011` | 0x1c3  |
    /// | 𝑋⁸ + 𝑋⁷ + 𝑋⁶ + 𝑋³ + 𝑋² + 𝑋  + 1 | `0b1_1100_1111` | 0x1cf  |
    /// | 𝑋⁸ + 𝑋⁷ + 𝑋⁶ + 𝑋⁴ + 𝑋² + 𝑋  + 1 | `0b1_1101_0111` | 0x1d7  |
    /// | 𝑋⁸ + 𝑋⁷ + 𝑋⁶ + 𝑋⁴ + 𝑋³ + 𝑋² + 1 | `0b1_1101_1101` | 0x1dd  |
    /// | 𝑋⁸ + 𝑋⁷ + 𝑋⁶ + 𝑋⁵ + 𝑋² + 𝑋  + 1 | `0b1_1110_0111` | 0x1e7  |
    /// | 𝑋⁸ + 𝑋⁷ + 𝑋⁶ + 𝑋⁵ + 𝑋⁴ + 𝑋  + 1 | `0b1_1111_0011` | 0x1f3  |
    /// | 𝑋⁸ + 𝑋⁷ + 𝑋⁶ + 𝑋⁵ + 𝑋⁴ + 𝑋² + 1 | `0b1_1111_0101` | 0x1f5  |
    /// | 𝑋⁸ + 𝑋⁷ + 𝑋⁶ + 𝑋⁵ + 𝑋⁴ + 𝑋³ + 1 | `0b1_1111_1001` | 0x1f9  |
    pub const IRREDUCTIBLE_POLYNOMIAL: u16 = 0x1cf;

    /// The null polynomial is neutral for addition
    pub const ZERO: Self = Self(0);

    /// Constant polynomial 1 is neutral for multiplication
    pub const ONE: Self = Self(1);

    /// Maps a byte to an element of this field
    ///
    /// # Parameters
    ///
    /// | Parameter | Description | Notes |
    /// | --------- | ----------- | ----- |
    /// | `value` | byte-mapped value of the element | |
    ///
    /// # Output
    ///
    /// Newly created element
    #[inline]
    pub const fn new(value: u8) -> Self {
        Self(value)
    }

    /// Retrieves the byte value of an element
    ///
    /// # Output
    ///
    /// Byte value of the element
    #[inline]
    pub const fn value(self) -> u8 {
        self.0
    }

    /// Quick reference table holding the multiplication table for all elements, calculated at compile time
    ///
    /// The nightly toolchain is required with `#![feature(const_eval_limit)]` and
    /// `#![const_eval_limit = "0"]` for this to work
    const ALL_PRODUCTS: [[Self; 256]; 256] = {
        let mut all_products = [[Self::new(0); 256]; 256];
        let mut i = 255_u8;
        while i > 0 {
            let mut j = 255_u8;
            while j > 0 {
                // Perform polynomial products.
                // Since all elements are polynomials, the product could exceed the maximal allowed degree (7).
                // The product is performed modulo the chosen irreductible polynomial
                let mut product: u8 = 0;
                let mut accumulator = i as u16;
                let mut rhs_remainder = j;

                while accumulator != 0 && rhs_remainder != 0 {
                    if rhs_remainder % 2 == 1 {
                        product ^= accumulator as u8;
                    }

                    rhs_remainder >>= 1;

                    accumulator <<= 1;
                    if accumulator >= 0b1_0000_0000 {
                        accumulator ^= Self::IRREDUCTIBLE_POLYNOMIAL;
                    }

                    if rhs_remainder == 0 {
                        break;
                    }
                }

                all_products[i as usize][j as usize] = Self::new(product);

                j -= 1;
            }
            i -= 1;
        }

        all_products
    };

    /// Quick reference table holding the multiplicative inverses of all elements, calculated at compile time
    ///
    /// - None for input 0
    /// - else an element that, multiplied with self, results in the neutral value for multiplication.
    ///
    /// The nightly toolchain is required with `#![feature(const_eval_limit)]` and
    /// `#![const_eval_limit = "0"]` for this to work
    const ALL_ELEMENT_INVERSES: [Option<Self>; 256] = {
        let mut all_inverses = [None; 256];
        let mut i = 255_u8;
        while i > 0 {
            let mut j: u8 = 255;
            while j > 0 {
                // cannot use == on Element, because this function is not const per Eq trait
                if Self::ALL_PRODUCTS[i as usize][j as usize].0 == Self::ONE.0 {
                    all_inverses[i as usize] = Some(Self::new(j));
                    break;
                }
                j -= 1;
            }
            i -= 1;
        }

        all_inverses
    };

    /// Quick reference table holding the powers table for all elements, calculated at compile
    /// time, up to the 256th power
    ///
    /// Per this table, 0⁰ is equal to 0.
    ///
    /// The nightly toolchain is required with `#![feature(const_eval_limit)]` and
    /// `#![const_eval_limit = "0"]` for this to work
    const ALL_POWERS: [[Self; 256]; 256] = {
        let mut all_powers = [[Self::new(0); 256]; 256];
        let mut i = 0_u8;
        loop {
            let mut accumulator = Self::ONE;
            let mut j = 0;
            loop {
                all_powers[i as usize][j as usize] = accumulator;
                if j == 255 {
                    break;
                }

                accumulator = Self::ALL_PRODUCTS[accumulator.value() as usize][i as usize];
                j += 1;
            }
            if i == 255 {
                break;
            }
            i += 1;
        }

        all_powers
    };

    /// Computes the multiplicative inverse of an element
    ///
    /// 0 has no multiplicative inverse, hence this function returns None for 0.
    ///
    /// # Output
    ///
    /// - None for input 0
    /// - else an element that, multiplied with self, results in the neutral value for multiplication.
    #[inline]
    pub const fn inverse(self) -> Option<Self> {
        match self {
            Self::ZERO => None,
            _ => Self::ALL_ELEMENT_INVERSES[self.value() as usize],
        }
    }

    /// Computes the exponentiation of an element by direct lookup into a table
    ///
    /// # Parameters
    ///
    /// | Parameter | Description | Notes |
    /// | --------- | ----------- | ----- |
    /// | `exponent` | exponent | |
    ///
    /// # Output
    ///
    /// - the neutral multiplication value for 0⁰
    /// - self to the power of `exponent`.
    #[inline]
    pub const fn power(self, exponent: u8) -> Self {
        Self::ALL_POWERS[self.value() as usize][exponent as usize]
    }
}

/// Computes a list of polynomials of degree 8, with coefficients in ℤ/2ℤ,
/// that cannot be expressed as the product of two lesser degree polynomials.
#[allow(dead_code)]
pub fn compute_irreductible_polynomials() -> Vec<u16> {
    let mut irreductible_polynomials = vec![];

    for poly in 0x100_u16..=0x1ffu16 {
        let mut factored = false;
        'poly: for i in 1..=255 {
            for j in i..=255 {
                let mut product = 0_u16;
                let mut accumulator = i as u16;
                let mut remainder = j as u16;
                while accumulator != 0 && accumulator < 0b10_0000_0000 && remainder != 0 {
                    if remainder % 2 == 1 {
                        product ^= accumulator;
                    }

                    accumulator <<= 1;
                    remainder >>= 1;
                }

                if product == poly {
                    factored = true;
                    break 'poly;
                }
            }
        }

        if !factored {
            irreductible_polynomials.push(poly);
        }
    }

    irreductible_polynomials
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn addition_is_commutative() {
        for i in 0..=255 {
            for j in 0..=255 {
                assert_eq!(
                    Element::new(i) + Element::new(j),
                    Element::new(j) + Element::new(i),
                    "{0} + {1} = {2:?}, {1} + {0} = {3:?}",
                    i,
                    j,
                    Element::new(i) + Element::new(j),
                    Element::new(j) + Element::new(i),
                );
            }
        }
    }

    #[test]
    fn addition_is_associative() {
        for i in 0..=255 {
            for j in 0..=255 {
                for k in 0..=255 {
                    let (i, j, k) = (Element::new(i), Element::new(j), Element::new(k));
                    assert_eq!(
                        i + (j + k),
                        (i + j) + k,
                        "{0:?} + ({1:?} + {2:?}) = {3:?}, ({0:?} + {1:?}) + {2:?} = {4:?}",
                        i,
                        j,
                        k,
                        i + (j + k),
                        (i + j) + k,
                    );
                }
            }
        }
    }

    #[test]
    fn zero_is_neutral_for_addition() {
        for i in 0..=255 {
            let i = Element::new(i);
            assert_eq!(
                i + Element::ZERO,
                i,
                "{0:?} + neutral = {1:?}",
                i,
                i + Element::ZERO,
            );
        }
    }

    #[test]
    fn addition_table_is_unique() {
        for i in 1..=255 {
            for j in i..=255 {
                for k in j..=255 {
                    if k == j {
                        continue;
                    }

                    let (i, j, k) = (Element::new(i), Element::new(j), Element::new(k));
                    assert_ne!(
                        i + j,
                        i + k,
                        "{0:?} + {1:?} = {0:?} + {2:?} = {3:?}",
                        i,
                        j,
                        k,
                        i + j,
                    );
                }
            }
        }
    }

    #[test]
    fn multiplication_is_commutative() {
        for i in 0..=255 {
            for j in 0..=255 {
                let (i, j) = (Element::new(i), Element::new(j));
                assert_eq!(
                    i * j,
                    j * i,
                    "{0:?} × {1:?} = {2:?}, {1:?} * {0:?} = {3:?}",
                    i,
                    j,
                    i * j,
                    j * i,
                );
            }
        }
    }

    #[test]
    fn multiplication_is_associative() {
        for i in 0..=255 {
            for j in 0..=255 {
                for k in 0..=255 {
                    let (i, j, k) = (Element::new(i), Element::new(j), Element::new(k));
                    assert_eq!(
                        i * (j * k),
                        (i * j) * k,
                        "{0:?} × ({1:?} × {2:?}) = {3:?}, ({0:?} × {1:?}) × {2:?} = {4:?}",
                        i,
                        j,
                        k,
                        i * (j * k),
                        (i * j) * k,
                    );
                }
            }
        }
    }

    #[test]
    fn one_is_neutral_for_multiplication() {
        for i in 0..=255 {
            let i = Element::new(i);
            assert_eq!(
                i * Element::ONE,
                i,
                "{0:?} × neutral = {1:?}",
                i,
                i * Element::ONE,
            );
        }
    }

    #[test]
    fn zero_is_sink_for_multiplication() {
        for i in 0..=255 {
            let i = Element::new(i);
            assert_eq!(
                i * Element::new(0),
                Element::new(0),
                "{0:?} × 0 = {1:?}",
                i,
                i * Element::new(0),
            );
        }
    }

    #[test]
    fn multiplication_table_is_unique() {
        for i in 1..=255 {
            for j in i..=255 {
                for k in j..=255 {
                    if k == j {
                        continue;
                    }
                    let (i, j, k) = (Element::new(i), Element::new(j), Element::new(k));
                    assert_ne!(
                        i * j,
                        i * k,
                        "{0:?} × {1:?} = {0:?} × {2:?} = {3:?}",
                        i,
                        j,
                        k,
                        i * j,
                    );
                }
            }
        }
    }

    #[test]
    fn subtraction_is_addition() {
        for i in 0..=255 {
            for j in 0..=255 {
                let (i, j) = (Element::new(i), Element::new(j));
                assert_eq!(
                    i + j,
                    i - j,
                    "{0:?} + {1:?} = {2:?}, {0:?} - {1:?} = {3:?}",
                    i,
                    j,
                    i + j,
                    i - j,
                );
            }
        }
    }

    #[test]
    fn element_is_self_opposite() {
        for i in 0..=255 {
            let i = Element::new(i);
            assert_eq!(-i, i, "-{0:?} = {1:?}", i, -i,);
        }
    }

    #[test]
    fn division_table_is_unique() {
        for i in 1..=255 {
            for j in i..=255 {
                for k in j..=255 {
                    if k == j {
                        continue;
                    }
                    let (i, j, k) = (Element::new(i), Element::new(j), Element::new(k));
                    assert_ne!(
                        i / j,
                        i / k,
                        "{0:?} ÷ {1:?} = {0:?} ÷ {2:?} = {3:?}",
                        i,
                        j,
                        k,
                        (i / j).unwrap_or_else(|| Element::new(0)),
                    );
                }
            }
        }
    }

    #[test]
    fn division_is_inverse_of_multiplication() {
        for i in 0..=255 {
            for j in 1..=255 {
                let (i, j) = (Element::new(i), Element::new(j));
                assert_eq!(
                    ((i * j) / j).unwrap(),
                    i,
                    "({0:?} × {1:?}) ÷ {1:?} ≠ {0:?}",
                    i,
                    j,
                );
                assert_eq!(
                    ((i / j).unwrap()) * j,
                    i,
                    "({0:?} ÷ {1:?}) × {1:?} ≠ {0:?}",
                    i,
                    j,
                );
            }
        }
    }

    #[test]
    fn division_is_last_power() {
        // There are 255 invertible elements. The power operation cycles every 255 operations, so
        // the 254th power must be the inverse.
        for i in 1..=255 {
            let i = Element::new(i);
            assert_eq!(
                i.power(254),
                i.inverse().unwrap(),
                "({0:?}^254 != {0:?}^-1)",
                i,
            );
        }
    }

    #[test]
    fn division_by_self_is_one() {
        for i in 1..=255 {
            let i = Element::new(i);
            assert_eq!(
                i / i,
                Some(Element::ONE),
                "{0:?} ÷ {1:?} = {2:?}",
                i,
                i,
                (i / i).unwrap(),
            );
        }
    }

    #[test]
    fn multiplication_is_distributive_over_addition() {
        for i in 0..=255 {
            for j in 0..=255 {
                for k in 0..=255 {
                    let (i, j, k) = (Element::new(i), Element::new(j), Element::new(k));
                    assert_eq!(
                        i * (j + k),
                        (i * j) + (i * k),
                        "{0:?} × ({1:?} + {2:?}) = {3:?}, {0:?} × {1:?} + {0:?} × {2:?} = {4:?}",
                        i,
                        j,
                        k,
                        i * (j + k),
                        (i * j) + (i * k),
                    );
                    assert_eq!(
                        (i + j) * k,
                        (i * k) + (j * k),
                        "({0:?} + {1:?}) × {2:?}) = {3:?}, {0:?} × {2:?} + {1:?} × {2:?} = {4:?}",
                        i,
                        j,
                        k,
                        (i + j) * k,
                        (i * k) + (j * k),
                    );
                }
            }
        }
    }

    #[test]
    fn power_recursion() {
        for i in 0..=255 {
            for k in 0..=254_u8 {
                let i = Element::new(i);
                assert_eq!(
                    i.power(k + 1),
                    i.power(k) * i,
                    "{0:?}^{1} × {0:?} = {3:?}, {0:?}^{2} = {4:?}",
                    i,
                    k,
                    k + 1,
                    i.power(k) * i,
                    i.power(k + 1)
                );
            }
        }
    }

    #[test]
    fn constant_irreductible_polynomial_is_indeed_irreductible() {
        for i in 1..=255 {
            for j in i..=255 {
                let mut product = 0_u16;
                let mut accumulator = i as u16;
                let mut remainder = j as u16;
                while accumulator != 0 && accumulator < 0b10_0000_0000 && remainder != 0 {
                    if remainder % 2 == 1 {
                        product ^= accumulator;
                    }

                    accumulator <<= 1;
                    remainder >>= 1;
                }

                assert_ne!(
                    product,
                    Element::IRREDUCTIBLE_POLYNOMIAL,
                    "{0:#b} × {1:#b} = {2:#b}",
                    i,
                    j,
                    Element::IRREDUCTIBLE_POLYNOMIAL
                );
            }
        }
    }

    #[bench]
    fn bench_addition(b: &mut test::Bencher) {
        b.iter(|| {
            for i in 0..=255 {
                for j in 0..=255 {
                    test::black_box(Element::new(i) + Element::new(j));
                }
            }
        });
    }

    #[bench]
    fn bench_opposite(b: &mut test::Bencher) {
        b.iter(|| {
            for i in 0..=255 {
                test::black_box(-Element::new(i));
            }
        });
    }

    #[bench]
    fn bench_subtraction(b: &mut test::Bencher) {
        b.iter(|| {
            for i in 0..=255 {
                for j in 0..=255 {
                    test::black_box(Element::new(i) - Element::new(j));
                }
            }
        });
    }

    #[bench]
    fn bench_multiplication(b: &mut test::Bencher) {
        b.iter(|| {
            for i in 0..=255 {
                for j in 0..=255 {
                    test::black_box(Element::new(i) * Element::new(j));
                }
            }
        });
    }

    #[bench]
    fn bench_inverse(b: &mut test::Bencher) {
        b.iter(|| {
            for i in 0..=255 {
                test::black_box(Element::new(i).inverse());
            }
        });
    }

    #[bench]
    fn bench_division(b: &mut test::Bencher) {
        b.iter(|| {
            for i in 0..=255 {
                for j in 0..=255 {
                    test::black_box(Element::new(i) / Element::new(j));
                }
            }
        });
    }

    #[bench]
    fn bench_power(b: &mut test::Bencher) {
        b.iter(|| {
            for i in 0..=255 {
                for j in 0..=255 {
                    test::black_box(Element::new(i).power(j));
                }
            }
        });
    }
}
