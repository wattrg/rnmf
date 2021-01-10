
/// ## Module containing mesh implementations
/// Each type of mesh is a sub-module within this module. The currently available mesh types are
/// * Cartesian
pub mod mesh;

/// ## Configuration module
/// Contains simulation wide settings which define the simulation
pub mod config;

/// ## Module defining the available boundary conditions
/// It defines the `BoundaryCondition` trait, which all DataFrames must have to fill the boundary
/// conditions. It also contains an enum `BCType` which lists the available boundary conditions.
pub mod boundary_conditions;

pub mod io;

#[cfg(feature="double_precision")]
pub type Real = f64;

#[cfg(not(feature="double_precision"))]
pub type Real = f32;
