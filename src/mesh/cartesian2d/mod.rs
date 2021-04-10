pub mod mesh;
pub mod dataframe;
pub mod block;

pub use crate::{Real, RealVec2, UIntVec2};

/// provides functionality to assist with indexing of data. This is separate to the index operator
trait Indexing{
    /// Returns the non-flat index of the item stored at a particular flattened location
    fn get_ij(&self, p: usize ) -> Option<(isize, isize, usize)>;

    /// Returns a borrow to the value stored at a particular non-flat index
    fn index(&self, i: isize, j: isize, n: usize) -> &Real;
}

/// function to convert a multi-dimensional coordinate into a single dimension coordinate
fn get_flat_index(i: isize, j: isize, n: usize, dimension: &UIntVec2, ng: usize) -> usize {

    let mut p = i + ng as isize;
    p += dimension[0] as isize * (j + ng as isize) 
            + (dimension[0] * dimension[1] * n) as isize;
    
    p as usize
}

/// Converts (2+1)D index into a 1D (flattened) index
fn get_ij_from_p(p: usize, dimension: &UIntVec2, n_nodes: usize, ng: usize) 
                                                        -> Option<(isize, isize, usize)>{
    if p >= n_nodes {
        Option::None
    }

    else{
        let i = (p % dimension[0]) as isize;
        let j;
        let d: isize;

        j = (p as isize - i)/dimension[0] as isize % dimension[1] as isize;

        d = (p as isize - j * dimension[1] as isize - i)
                    /dimension[0] as isize % 2;
        Some((i - ng as isize, j - ng as isize, d as usize))
    }
}

