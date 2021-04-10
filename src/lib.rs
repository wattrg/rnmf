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

// /// Module for handling global values
// pub mod global;

// /// Module for handling the execution
// pub mod solver;

/// Alias for either f32/f64 depending on if `disable_double` feature is enabled
#[cfg(not(feature="double_precision"))]
pub type Real = f64;

/// Alias for either f32/f64 depending on if `disable_double` feature is enabled
#[cfg(feature="disable_double")]
#[derive(Default)]
pub type Real = f32;



/// Array containing two real numbers
#[derive(Debug, Clone, Default)]
pub struct RealVec2(pub [Real; 2]);
impl std::ops::Index<usize> for RealVec2{
    type Output = Real;
    fn index(&self, indx: usize)->&Self::Output{
        &self.0[indx]
    }
}
impl std::ops::IndexMut<usize> for RealVec2{
    fn index_mut(&mut self, indx: usize) -> &mut Self::Output {
        &mut self.0[indx]
    }
}
impl std::ops::Deref for RealVec2{
    type Target = [Real; 2];

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

/// Array containing three real numbers
#[derive(Debug, Clone, Default)]
pub struct RealVec3(pub [Real; 3]);
impl std::ops::Index<usize> for RealVec3{
    type Output = Real;
    fn index(&self, indx: usize)->&Self::Output{
        &self.0[indx]
    }
}
impl std::ops::IndexMut<usize> for RealVec3{
    fn index_mut(&mut self, indx: usize) -> &mut Self::Output {
        &mut self.0[indx]
    }
}
impl std::ops::Deref for RealVec3{
    type Target = [Real; 3];

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}


/// Array containing two unsigned integers
#[derive(Debug, Clone, Default)]
pub struct UIntVec2(pub [usize; 2]);
impl std::ops::Deref for UIntVec2{
    type Target = [usize; 2];

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}
impl std::iter::FromIterator<usize> for UIntVec2{
    fn from_iter<I: IntoIterator<Item=usize>>(iter: I) -> Self {
        let mut result = UIntVec2([0, 0]);

        for (i, value) in iter.into_iter().enumerate(){
            result[i] = value;
        } 
        result
    }
}
impl std::ops::Index<usize> for UIntVec2{
    type Output = usize;
    fn index(&self, indx: usize)->&Self::Output{
        &self.0[indx]
    }
}
impl std::ops::IndexMut<usize> for UIntVec2{
    fn index_mut(&mut self, indx: usize) -> &mut Self::Output {
        &mut self.0[indx]
    }
}



/// Array containing three unsigned integers
#[derive(Debug, Clone, Default)]
pub struct UIntVec3(pub [usize; 3]);
impl std::ops::Deref for UIntVec3{
    type Target = [usize; 3];

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}
impl std::iter::FromIterator<usize> for UIntVec3{
    fn from_iter<I: IntoIterator<Item=usize>>(iter: I) -> Self {
        let mut result = UIntVec3([0, 0, 0]);

        for (i, value) in iter.into_iter().enumerate(){
            result[i] = value;
        } 
        result
    }
}
impl std::ops::Index<usize> for UIntVec3{
    type Output = usize;
    fn index(&self, indx: usize)->&Self::Output{
        &self.0[indx]
    }
}
impl std::ops::IndexMut<usize> for UIntVec3{
    fn index_mut(&mut self, indx: usize) -> &mut Self::Output {
        &mut self.0[indx]
    }
}




/// Slice containing two integers
#[derive(Debug, Clone, Default)]
pub struct IntVec2(pub [isize; 2]);
impl std::ops::Deref for IntVec2{
    type Target = [isize; 2];

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}
impl std::ops::Index<usize> for IntVec2{
    type Output = isize;
    fn index(&self, indx: usize)->&Self::Output{
        &self.0[indx]
    }
}
impl std::ops::IndexMut<usize> for IntVec2{
    fn index_mut(&mut self, indx: usize) -> &mut Self::Output {
        &mut self.0[indx]
    }
}

/// Slice containing three integers
#[derive(Debug, Clone, Default)]
pub struct IntVec3(pub [isize; 3]);
impl std::ops::Deref for IntVec3{
    type Target = [isize; 3];

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}
impl std::ops::Index<usize> for IntVec3{
    type Output = isize;
    fn index(&self, indx: usize)->&Self::Output{
        &self.0[indx]
    }
}
impl std::ops::IndexMut<usize> for IntVec3{
    fn index_mut(&mut self, indx: usize) -> &mut Self::Output {
        &mut self.0[indx]
    }
}




