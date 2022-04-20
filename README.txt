(C) 2022 Benoît Mortgat

shamir_secret_sharing - Rust implementation of Shamir secret sharing for use
with files.

This program has been tested on GNU/Linux and Windows.

USAGE
-----

Creation of Shamir secret shares:
    shamir_secret_sharing share <secret_file> <share_count> <threshold>

Where:
    * <secret_file> : sensitive file to be converted into Shamir shares
    * <share_count> : number of shares to create (2 to 255)
    * <threshold>   : minimal number of shares that are necessary in order
                      to rebuild the secret file.

Recovery of secret from shares:
    shamir_secret_sharing recover <recovered_file> <share1> <share2> ...

EXAMPLE
-------

You need a nightly rust toolchain to build. This crate makes heavy use of
constant lookup tables that are built at compile time.

In order to use on this file (example command-lines using a bash shell):

This example will make 5 shares with a threshold of 3 shares for recovery.

  # build the project
  $ cargo +nightly build --release

  # make the shares
  $ cargo +nightly run --release --  share   README.txt 5 3

  # rebuild the original file using only 3 shares picked at random
  $ cargo +nightly run --release --  recover README.txt{.out,_001,_004,_005}

  # compare original and recovered version: identical
  $ diff README.txt{,.out}

INSTALLATION
------------

You can use:
  $ cargo +nightly install                                        \
          --git https://github.com/salsifis/shamir_secret_sharing \
          --branch main                                           \
          shamir_secret_sharing

FAQ
---

Q. What is the share file format?
A. The file format consists of:
   * A 38-byte header:
     * 16 bytes with a fixed UUID that identifies the file format
     * 4 bytes with a file format version (currently 1)
     * 16 bytes with an UUID that identifies a set of shares
     * 1 byte storing the value of <threshold>
     * 1 byte storing a value specific to the share.
   * Content bytes (same byte length as the secret file).

Q. What are the mathematical objects used here?
A. The Galois field used is GF(256) using x^8 + x^7 + x^6 + x^3 + x^2 + x + 1
   as the irreductible polynomial for multiplication purposes. A list of
   suitable irreductible polynomials can be found in the source code.

Q. Any security weaknesses?
A. - Operations in the Galois field use lookup tables. As a consequence, they
     may be not run in constant time. When generating secret shares, be sure to
     run this program in a controlled environment that cannot be subject to
     timing attacks.
   - This uses thread_rng() for generation of random polynomials. If somehow
     thread_rng() non-predictability assumption was somehow broken in the
     future, then a single share may suffice to recover the secret.
   - There could be weaknesses in the way you use Shamir Secret Sharing.
     You should regularly audit your protocols for sharing and recovering
     sensitive data.

RUST
----

The source passes cargo clippy, and is auto-formatted with cargo fmt.

The sources include:
 * Tests (for use with cargo test)
 * Benchmarks (for use with cargo bench)
 * Rust documentation (for use with cargo doc or cargo rustdoc)
