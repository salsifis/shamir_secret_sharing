//! This program performs Shamir Secret Sharing, one file at a time
//!
//! Secret sharing works by converting a **secret** value to 𝑁 **shares**.
//!
//! Shares are individually useless, but any combination of ⩾ 𝑇 shares will
//! give the secret back, where 𝑇 is a value between 2 and 𝑁.
//!
//! Any combination of less than 𝑇 shares provides no indication at all regarding the secret
//! (all possible secret values are equiprobable).
//!
//! This can be used in order to protect a sensitive asset, such as a critical private key,
//! by “sharing” it among 10 individuals, and requiring that at least 5 of them combine their
//! shares together in order to be able to use the asset.
//!
//! In this implementation, the secret is a file, and each share is a file, that has:
//! - A fixed size header of 38 bytes
//!   - 16 bytes for a UUID (identifying this program)
//!   - 4 bytes for a file format version, currently 1
//!   - 16 bytes for a UUID, identifying the set of shares
//!   - 1 byte for 𝑇
//!   - 1 byte specific to the share, not equal to 0
//! - Contents with the same byte size as the secret
//!
//! This program can perform two operations:
//!
//! - Read an input, secret file, and create 𝑁 shares given a parameter 𝑇;
//! - Read at least 𝑇 shares, and output the recovered secret file.
//!

#![feature(const_eval_limit)]
#![feature(test)]
#![const_eval_limit = "0"]
// enabled clippy linters categories
#![warn(
    clippy::all,
    clippy::pedantic,
    clippy::nursery,
    clippy::cargo,
    clippy::style
)]
// disabled individual clippy linters
#![allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]

/// Dilutes secrets in a secure way
///
/// This program performs Shamir's secret sharing on files.
///
/// Shamir's secret sharing consists of converting a sensitive file (secret file)
/// into a number N of shares (also files).
///
/// The only thing an individual share reveals about the secret is its file size.
///
/// It is possible to combine T shares (or more), where T can be equal to N or less than N,
/// in order to be able to rebuild the original secret file.
///
/// Any combination of less than T shares provides no indication at all regarding the secret
/// (all possible secret files of the same size are equiprobable).
///
/// This can be used in order to protect a sensitive asset, such as a critical private key,
/// by "sharing" it among 10 individuals, and requiring that at least 5 of them combine their
/// shares together in order to be able to use the asset.
#[derive(StructOpt)]
enum ShamirMainCommand {
    #[structopt(name = "share")]
    /// Transform a secret file into shares
    Share(ShareCommand),

    #[structopt(name = "recover")]
    /// Rebuild a secret file from shares
    Recover(RecoverCommand),
}

/// Converts a sensitive file (secret file) into a number N of shares (also files).
///
/// A parameter named T indicates how many shares must be put in common in order to
/// rebuild the original file. T must be more than 1, and less than or equal to N.
///
/// The output files will have the same path as the input plus the suffixes _001 to _255.
#[derive(StructOpt)]
struct ShareCommand {
    /// Path to the input secret file
    secret_file: String,
    /// Number of output shares (value of N)
    num_shares: u8,
    /// Number of shares to be used for recovery of secret (value of T)
    threshold: u8,
}

/// Rebuilds a sensitive file given a list of share files.
///
/// This operation could fail if there are not enough share files, or if
/// a share has been used more than once.
#[derive(StructOpt)]
struct RecoverCommand {
    /// Path to the output secret file
    secret_file: String,
    /// Paths to the input share files
    shares: Vec<String>,
}

mod galois_field_256;
mod polynomials;
mod shamir;

/// A UUID inserted in the header of each share, in binary form.
///
/// This works as a magic value, identifying that a share was really created by this program
const PROGRAM_UUID: uuid::Uuid =
    uuid::Uuid::from_u128(0x1bc0_96ca_eec2_4565_8634_f34a_5903_fcda_u128);

/// A version number for the format of a share.
const FILE_FORMAT_VERSION: u32 = 1;

use structopt::StructOpt;

/// Reads a secret file, and outputs 𝑁 shares.
///
/// Here we define:
///
/// * 𝑁: the number of shares created from the secret:
///   * 𝑁 must be between 2 and 255, as a single share would reveal the secret.
/// * 𝑇: the minimal number of shares (threshold) that must be combined to recover the secret:
///   * 𝑇 must be less than or equal to 𝑁
///   * 𝑇 must be greater than or equal to 2, as a threshold of 1 would reveal the secret.
///
/// # Parameters
///
/// | Parameter | Description | Notes |
/// | --------- | ----------- | ----- |
/// | `filepath` | Path to the secret file | |
/// | `num_shares` | Number of output files to create (𝑁) | Must be between 2 and 255 |
/// | `threshold` | Minimal number of output files to combine for the secret to be recovered (𝑇) | Must be between 2 and 𝑁 |
///
/// # File system
///
/// Creates 𝑁 files with the same path as `filepath` and suffixes `_001`..`_𝑁`.
///
/// # Output
///
/// - Ok without error
/// - Err(...) otherwise
fn split_secret_file_into_shares(
    filepath: &str,
    num_shares: u8,
    threshold: u8,
) -> Result<(), Box<dyn std::error::Error>> {
    use std::fs::File;
    use std::io::prelude::*;

    let unique_uuid_for_this_operation = uuid::Uuid::new_v4();

    if num_shares < 2 {
        return Err("Must have at least two shares".into());
    }
    if threshold < 2 {
        return Err("Must have a threshold of at least two".into());
    }
    if threshold > num_shares {
        return Err("Must have at least as many shares as the threshold for recovery".into());
    }

    let mut rng = rand::thread_rng();
    let abscissae = shamir::create_distinct_x_values(num_shares, &mut rng);

    let input_file = std::io::BufReader::new(File::open(filepath)?);
    let mut output_files: Vec<_> = vec![];
    for index in 1..=num_shares {
        let file_name = format!("{}_{:03}", filepath, index);
        let opened_file = std::io::BufWriter::new(File::create(file_name)?);
        output_files.push(opened_file);
    }

    if abscissae.len() != output_files.len() {
        return Err("Logic error".into());
    }

    // File header
    // - Guid for the program
    // - Version of the share file format
    // - Guid for this secret sharing (prevents mixing different secrets together)
    // - 𝑇
    // - 𝑁
    // - an abscissa used for all polynomial interpolations
    for (idx, output) in output_files.iter().enumerate() {
        bincode::serialize_into(
            output.get_ref(),
            &(
                PROGRAM_UUID.as_bytes(),
                FILE_FORMAT_VERSION,
                unique_uuid_for_this_operation.as_bytes(),
                threshold,
                abscissae[idx],
            ),
        )?;
    }

    for input_byte in input_file.bytes() {
        let shares =
            shamir::create_secret_shares_from_value(input_byte?, &abscissae, threshold, &mut rng);
        if shares.len() != output_files.len() {
            return Err("Logic error".into());
        }
        for index in 0..shares.len() {
            output_files[index].write_all(&[shares[index]])?;
        }
    }
    Ok(())
}

/// Reads at least 𝑇 shares, and outputs the secret file
///
/// # Parameters
///
/// | Parameter | Description | Notes |
/// | --------- | ----------- | ----- |
/// | `inputs` | Paths to the shares | |
/// | `output` | Path to the recovered file | |
///
/// # File system
///
/// Creates `output`
///
/// # Output
///
/// - Ok without error
/// - Err(...) otherwise
fn recover_shared_secret(
    inputs: Vec<String>,
    output: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    use std::fs::File;
    use std::io::prelude::*;

    let mut input_files: Vec<_> = Vec::new();
    for input_file_name in inputs {
        let opened_file = std::io::BufReader::new(File::open(input_file_name)?);
        input_files.push(opened_file);
    }

    let mut output_file = std::io::BufWriter::new(File::create(output)?);

    let mut program_uuids: Vec<uuid::Uuid> = Vec::new();
    for stream in &input_files {
        program_uuids.push(uuid::Uuid::from_bytes(bincode::deserialize_from(
            stream.get_ref(),
        )?));
    }
    let program_uuids = program_uuids;
    if program_uuids.into_iter().any(|uuid| uuid != PROGRAM_UUID) {
        return Err("One of the input files was not created by this program".into());
    }

    let mut file_format_versions: Vec<u32> = Vec::new();
    for stream in &input_files {
        file_format_versions.push(bincode::deserialize_from(stream.get_ref())?);
    }
    let file_format_versions = file_format_versions;
    if file_format_versions
        .into_iter()
        .any(|ver| ver != FILE_FORMAT_VERSION)
    {
        return Err(
            "One of the input files was created by a non-supported version of this program".into(),
        );
    }

    let mut file_identifiers: Vec<uuid::Uuid> = Vec::new();
    for stream in &input_files {
        file_identifiers.push(uuid::Uuid::from_bytes(bincode::deserialize_from(
            stream.get_ref(),
        )?));
    }
    file_identifiers.sort();
    file_identifiers.dedup();
    if file_identifiers.len() != 1 {
        return Err("The secrets shares do not pertain to the same set of secrets".into());
    }

    let mut threshold: Vec<u8> = Vec::new();
    for stream in &input_files {
        threshold.push(bincode::deserialize_from(stream.get_ref())?);
    }
    threshold.sort_unstable();
    threshold.dedup();
    if threshold.len() != 1 {
        return Err("The secrets shares do not pertain to the same set of secrets".into());
    }
    let threshold = threshold[0];

    let mut abscissae: Vec<u8> = Vec::new();
    for stream in &input_files {
        abscissae.push(bincode::deserialize_from(stream.get_ref())?);
    }
    let abscissae = abscissae;
    let mut sorted_abscissae = abscissae.clone();
    sorted_abscissae.sort_unstable();
    sorted_abscissae.dedup();
    if sorted_abscissae.len() != input_files.len() {
        return Err("The same share has been given twice".into());
    }

    let mut y_values_as_options: Vec<Option<u8>> = Vec::with_capacity(input_files.len());
    let mut y_values: Vec<u8> = Vec::with_capacity(input_files.len());
    loop {
        y_values_as_options.clear();
        y_values_as_options.extend(input_files.iter_mut().map(|stream| {
            let mut x = [0];
            match stream.read(&mut x) {
                Ok(1) => Some(x[0]),
                _ => None,
            }
        }));

        if y_values_as_options.iter().any(|x| *x == None) {
            if y_values_as_options.iter().all(|x| *x == None) {
                break;
            }
            return Err("An input file has reached its end".into());
        }

        y_values.clear();
        y_values.extend(y_values_as_options.iter().map(|x| match x {
            Some(val) => *val,
            _ => 0,
        }));

        let recovered_byte = shamir::recover_secret_from_points(&abscissae, &y_values, threshold);
        output_file.write_all(&[recovered_byte])?;
    }

    Ok(())
}

/// Main program
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let operation = ShamirMainCommand::from_args();
    match operation {
        ShamirMainCommand::Share(cmd) => {
            split_secret_file_into_shares(&cmd.secret_file, cmd.num_shares, cmd.threshold)?;
        }
        ShamirMainCommand::Recover(cmd) => recover_shared_secret(cmd.shares, &cmd.secret_file)?,
    };
    Ok(())
}
