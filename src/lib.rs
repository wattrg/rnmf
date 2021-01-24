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


#[derive(Debug, Clone)]
pub struct RealVec1(pub [Real; 1]);

/// Array containing two real numbers
#[derive(Debug, Clone)]
pub struct RealVec2(pub [Real; 2]);
impl std::ops::Index<usize> for RealVec2{
    type Output = Real;
    fn index(&self, indx: usize)->&Self::Output{
        &self.0[indx]
    }
}

/// Array containing three real numbers
#[derive(Debug, Clone)]
pub struct RealVec3(pub [Real; 3]);

#[derive(Debug, Clone)]
pub struct UIntVec1(pub [usize; 1]);

/// Array containing two unsigned integers
#[derive(Debug, Clone)]
pub struct UIntVec2(pub [usize; 2]);
impl std::ops::Index<usize> for UIntVec2{
    type Output = usize;
    fn index(&self, indx: usize)->&Self::Output{
        &self.0[indx]
    }
}

/// Array containing three unsigned integers
#[derive(Debug, Clone)]
pub struct UIntVec3(pub [usize; 3]);

#[derive(Debug, Clone)]
pub struct IntVec1(pub [isize; 1]);

/// Slice containing two integers
#[derive(Debug, Clone)]
pub struct IntVec2(pub [isize; 2]);
impl std::ops::Index<usize> for IntVec2{
    type Output = isize;
    fn index(&self, indx: usize)->&Self::Output{
        &self.0[indx]
    }
}

/// Slice containing three integers
#[derive(Debug, Clone)]
pub struct IntVec3(pub [isize; 3]);

/// Contains various types of data which can be used
#[derive(Debug, Clone)]
pub enum RnmfType{
    Usize(usize),
    Isize(isize),
    Real(Real),
    RealVec1(RealVec1),
    RealVec2(RealVec2),
    RealVec3(RealVec3),
    UIntVec2(UIntVec2),
    UIntVec3(UIntVec3),
    IntVec2(IntVec2),
    IntVec3(IntVec3),
}

