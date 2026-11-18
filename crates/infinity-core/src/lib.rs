//! Infinity-Win: progressive, cluster-first, position-aware comparison of nearly
//! identical genomes.
//!
//! # Layout
//!
//! The crate is a library (section 0 rule 2: the engine must never depend on
//! the GUI). Everything a consumer can see arrives as an [`events::Event`];
//! the command-line tool and the desktop program are two consumers of the same
//! stream.
//!
//! | Module | Section of the outline |
//! |---|---|
//! | [`params`] | 6.1 parameters and the experiment switches |
//! | [`seq`] | 6.2 encoding, the anchor hash, level thresholds |
//! | [`fasta`] | 6.2.1 input, 9 accepted files |
//! | [`genome`] | 6.2 scan, per-genome statistics |
//! | [`blocks`] | 6.3 nested blocks and their strand-independent hashes |
//! | [`sketch`] | 6.4 sketches and the sketch distance |
//! | [`cluster`] | 6.4 greedy clustering and the backbone tree |
//! | [`tree`] | 6.4/6.9 neighbour joining, midpoint rooting, Newick, Robinson-Foulds |
//! | [`pivot`] | 6.5 the pivot index |
//! | [`compare`] | 6.6 comparing a member to its pivot |
//! | [`align`] | 6.6 D and 6.7 banded affine alignment, x-drop extension |
//! | [`variant`] | 6.7 resolving a difference piece |
//! | [`coverage`] | 6.8 what a member does not cover |
//! | [`profile`] | 6.9 profiles, pairwise distances, cluster trees |
//! | [`engine`] | 6.11 the order of work |
//! | [`events`] | 7 the message stream |
//! | [`output`] | 9 the files in the output folder |
//!
//! # Coordinates
//!
//! 0-based and half-open `[start, end)` inside the program. Output files use
//! 1-based positions, converted once in [`output`].
//!
//! # Determinism
//!
//! Same inputs and parameters give byte-identical output files whatever the
//! thread count (section 0 rule 3). That constrains the whole design: results
//! are computed in parallel but always merged in a fixed order, no iteration
//! over an unordered map ever reaches the output, and every tie is broken by a
//! documented rule.

#![forbid(unsafe_code)]
#![warn(missing_debug_implementations)]

pub mod align;
pub mod backbone;
pub mod blocks;
pub mod cluster;
pub mod compare;
pub mod counters;
pub mod engine;
pub mod events;
pub mod export;
pub mod fasta;
pub mod filter;
pub mod genome;
pub mod interval;
pub mod output;
pub mod params;
pub mod pivot;
pub mod profile;
pub mod seq;
pub mod sketch;
pub mod summary;
pub mod tree;
pub mod variant;

pub use params::Parameters;

/// Version reported by `infinity --version` and written into `run.json`.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

/// Name of the program, used in output files and messages.
pub const PROGRAM: &str = "infinity";
