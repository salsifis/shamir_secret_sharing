//! This module performs Shamir Secret Sharing, one byte at a time
//!
//! Secret sharing works by converting a **secret** value to 𝑁 **shares**.
//!
//! Shares are individually useless, but any combination of ⩾ 𝑇 shares will
//! give the secret back.
//!
//! Any combination of < 𝑇 shares provides no indication at all regarding the secret
//! (all possible secret values are equiprobable).
//!
//! In this module we define:
//!
//! * 𝑠: the secret value
//! * 𝑁: the number of shares created from the secret:
//!   * 𝑁 must be between 2 and 255, as a single share would reveal the secret.
//! * 𝑇: the minimal number of shares (threshold) that must be combined to recover the secret:
//!   * 𝑇 must be less than or equal to 𝑁
//!   * 𝑇 must be greater than or equal to 2, as a threshold of 1 would reveal the secret.
//!
use crate::polynomials::GF256Polynomial;

pub type Byte = u8;

/// Creates a set of 𝑁 distinct, random abscissae for creation of Shamir shares
///
/// This is a prerequisite for creation of a Shamir share.
///
/// Since this module works byte by byte, it is possible to reuse the abscissae that this
/// function outputs for another byte as long as the same polynomial is not re-used.
///
/// # Parameters
///
/// | Parameter | Description | Notes |
/// | --------- | ----------- | ----- |
/// | `number_of_abscissae` | Desired cardinality of output (𝑁) | Must be between 2 and 255 |
/// | `rng` | A random number generator | |
///
/// # Output
///
/// Vector of 𝑁 distinct values. All values are in the range 1..=255 inclusive
pub fn create_distinct_x_values<RNG>(number_of_abscissae: u8, rng: &mut RNG) -> Vec<Byte>
where
    RNG: rand::Rng,
{
    use rand::prelude::IteratorRandom;
    (1_u8..=255_u8).choose_multiple(&mut *rng, number_of_abscissae as usize)
}

/// Creates 𝑁 Shamir shares with recovery threshold 𝑇 for a secret byte 𝑠.
///
/// Shamir shares are points (𝑥, 𝑦). The values of 𝑥 are provided to this function.
///
/// Any subset of distinct share, with a cardinality of at minimum 𝑇, can be used to recover
/// the secret byte.
///
/// This works by creating a polynomial of degree 𝑇-1, such the secret value is the constant
/// coefficient (of degree zero) of the polynomial. This means the point (0, 𝑠) is along
/// the curve of the polynomial.
///
/// The input of this function is (𝑥ᵢ), the output is (𝑦ᵢ) such that all points (𝑥ᵢ, 𝑦ᵢ) are
/// along the curve of the polynomial.
///
/// The recovery of the secret necessitates at least 𝑇 distinct points along its curve, and
/// the knowledge of the degree of the polynomial.
///
/// Using this interpolation, the value when 𝑥 = 0 provides the secret.
///
/// # Parameters
///
/// | Parameter | Description | Notes |
/// | --------- | ----------- | ----- |
/// | `secret_value` | The secret byte 𝑠 | |
/// | `abscissae_of_secret_shares` | 𝑁 distinct values between 1 and 255 inclusive | Use `create_distinct_x_values` to create these |
/// | `min_shares_for_recovery` | value of 𝑇 | This must be at least 2, otherwise the polynomial would be constant and 𝑠 would be in plain text |
/// | `rng` | A random number generator | This is used for creating the random polynomial
///
/// # Output
///
/// Vector of 𝑁 distinct ordinates that can be zipped with `abscissae_of_secret_shares` to represent points along a random polynomial.
///
/// Reconstitution of the secret needs at least 𝑇 points (abscissa, ordinate).
pub fn create_secret_shares_from_value<RNG>(
    secret_value: Byte,
    abscissae_of_secret_shares: &[Byte],
    min_shares_for_recovery: u8,
    rng: &mut RNG,
) -> Vec<Byte>
where
    RNG: rand::Rng,
{
    let number_of_shares = abscissae_of_secret_shares.len();

    assert!(number_of_shares <= 255, "Too many shares");

    let number_of_shares = number_of_shares as u8;

    assert!(
        number_of_shares >= 2,
        "At least two shares are necessary for Shamir secret sharing"
    );
    assert!(
        min_shares_for_recovery <= number_of_shares,
        "There must be more shares than the threshold for recovery"
    );
    assert!(
        min_shares_for_recovery >= 2,
        "Using less than 2 shares for recovery is useless"
    );
    assert!(
        abscissae_of_secret_shares.iter().all(|x| *x != 0),
        "The abscissa 0 is the plain secret"
    );

    let degree = min_shares_for_recovery - 1;

    let random_polynomial = GF256Polynomial::new_random(degree, Some(secret_value), &mut *rng);

    abscissae_of_secret_shares
        .iter()
        .map(|x| random_polynomial.value_at(*x))
        .collect()
}

/// Recover a shamir secret from at least 𝑇 shares
///
/// # Parameters
///
/// | Parameter | Description | Notes |
/// | --------- | ----------- | ----- |
/// | `x_values` | at least 𝑇 values between 1 and 255 inclusive | |
/// | `y_values` | at least 𝑇 values | As many values as `x_values` |
/// | `min_shares_for_recovery` | actual value of 𝑇 | Since `x_values` and `y_values` can hold more than 𝑇 points along the polynomial, this indicates how many points will be really used |
///
/// # Output
///
/// value of 𝑠.
pub fn recover_secret_from_points(
    x_values: &[Byte],
    y_values: &[Byte],
    min_shares_for_recovery: u8,
) -> Byte {
    assert!(
        min_shares_for_recovery >= 2,
        "Using less than 2 shares for recovery is useless"
    );
    assert!(
        x_values.len() == y_values.len(),
        "Should have as many abscissae as ordinates"
    );
    assert!(
        x_values.len() >= min_shares_for_recovery.into(),
        "Not enough points provided"
    );

    GF256Polynomial::new_from_points(x_values, y_values, min_shares_for_recovery - 1).value_at(0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn shamir_secret_is_recoverable() {
        // traits providing sample_iter
        use rand::distributions::Distribution;

        let mut rng = rand::thread_rng();

        // Take 5 distinct num_shares values
        let tested_num_shares = rand::distributions::Uniform::new_inclusive(2_u8, 255_u8)
            .sample_iter(&mut rng)
            .take(5)
            .collect::<Vec<u8>>();

        for num_shares in tested_num_shares {
            // Take 5 threshold values
            let tested_thresholds = rand::distributions::Uniform::new_inclusive(2_u8, num_shares)
                .sample_iter(&mut rng)
                .take(2)
                .collect::<Vec<u8>>();

            for threshold in tested_thresholds {
                // take 5 secrets
                let attempted_secrets = rand::distributions::Uniform::new_inclusive(0_u8, 255_u8)
                    .sample_iter(&mut rng)
                    .take(5)
                    .collect::<Vec<u8>>();

                for secret in attempted_secrets {
                    let x_values = create_distinct_x_values(num_shares, &mut rng);

                    let points =
                        create_secret_shares_from_value(secret, &x_values, threshold, &mut rng);

                    assert_eq!(points.len(), num_shares.into());

                    let recovered_secret =
                        recover_secret_from_points(&x_values, &points, threshold);

                    assert_eq!(recovered_secret, secret);
                }
            }
        }
    }
}
