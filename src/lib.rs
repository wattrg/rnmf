//! `rnmf` aims to provide functionality to assist in writing efficient numerical codes for a
//! variety of applications.
//! 
//! # Example
//! I will put a simple example here



/// ## Module containing mesh implementations
/// Each type of mesh is a sub-module within this module. The currently available mesh types are
/// * Cartesian
pub mod mesh;

/// Contains simulation wide settings which define the simulation
pub mod config;

/// Module defining available boundary conditions
pub mod boundary_conditions;

/// Module to assist with input and output 
pub mod io;

/// Alias for either f32/f64 depending on if `disable_double` feature is enabled
#[cfg(not(feature="double_precision"))]
pub type Real = f64;

/// Alias for either f32/f64 depending on if `disable_double` feature is enabled
#[cfg(feature="disable_double")]
pub type Real = f32;

/// Slice containing two real numbers
pub type RealVec2 = [Real; 2];

/// Slice containing three real numbers
pub type RealVec3 = [Real; 3];

/// Slice containing two unsigned integers
pub type UIntVec2  = [usize; 2];

/// Slice containing three unsigned integers
pub type UIntVec3  = [usize; 3];

/// Slice containing two integers
pub type IntVec2  = [isize; 2];

/// Slice containing three integers
pub type IntVec3  = [isize; 3];
