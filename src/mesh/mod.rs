use core::ops::Index;
use crate::Real;

//pub trait State<S=Self>
//where State<S>:
//      std::ops::Add
//    + std::ops::Sub
//    + std::ops::Mul
//{
//
//}

//impl State for Real{}

/// # Two dimensional cartesian mesh
///
/// A mesh where each cell is the same size, and has edges perpendicular to the *x* and *y* axis.
pub mod cartesian2d;

/// Allows for operations to be performed on any block, regardless of the type of 
/// mesh contained within
pub trait Block: Index<String>{}

/// Allows for operations to be performed on any dataframe, regardless of the type of 
/// mesh stored within
pub trait DataFrame{}


/// # Data structure storing all mesh information
/// `Domain` is the data structure describing the domain, and all the information stored inside.
/// The structure itself is a light wrapper around std::Vec, but with some additional methods
/// defined to assist with creating domains, and 
/// automatically performing performing some operations.
/// 
/// The domain consists of one or more `Blocks`, each
/// containing a portion of the mesh. The domain may be split into separate `Blocks` to represent
/// complex geometry, or to distribute work among multiple threads. Each `Block` has at most 
/// one thread assigned to it, but each thread may be assigned to multiple `Block`s.
#[derive(Debug)]
pub struct Domain<B: Block> {
    blocks: Vec<B>,
}

