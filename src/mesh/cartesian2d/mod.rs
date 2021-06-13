pub use crate::{Real, RealVec2, UIntVec2};
use super::Domain;
use crate::boundary_conditions::BCs;
use block::CartesianBlock;

/// 
pub mod mesh;
pub mod dataframe;
pub mod block;


/// provides functionality to assist with indexing of data. This is separate to the index operator
trait Indexing{
    /// Returns the non-flat index of the item stored at a particular flattened location
    fn get_ij(&self, p: usize ) -> Option<(isize, isize, usize)>;

    /// Returns a borrow to the value stored at a particular non-flat index
    fn index(&self, i: isize, j: isize, n: usize) -> &Real;
}

/// function to convert a multi-dimensional coordinate into a single dimension coordinate
fn get_flat_index(i: isize, j: isize, n: usize, dimension: &UIntVec2, ng: usize) -> Option<usize> {
    // make sure i, j is a valid point in the mesh
    if (i < -(ng as isize)) | (i >= (dimension[0] + ng) as isize) 
     | (j < -(ng as isize)) | (j >= (dimension[1] + ng) as isize) {
        return None;
    }
    Some(get_flat_index_unchecked(i,j,n,dimension,ng))
}

fn get_flat_index_unchecked(i: isize, j: isize, n:usize, dimension: &UIntVec2, ng: usize) -> usize{
    let mut p = i + ng as isize;
    p += dimension[0] as isize * (j + ng as isize) 
            + (dimension[0] * dimension[1] * n) as isize;
    
    p as usize
}

/// Converts (2+1)D index into a 1D (flattened) index
fn get_ij_from_p(p: usize, dimension: &UIntVec2, n_nodes: usize, ng: usize) 
                                         -> Option<(isize, isize, usize)>{
    if p > n_nodes {
        return Option::None
    }

    let i = (p % dimension[0]) as isize;
    let j;
    let d: isize;

    j = (p as isize - i)/dimension[0] as isize % dimension[1] as isize;

    d = (p as isize - j * dimension[1] as isize - i) / (dimension[0]*dimension[1]) as isize % 2;
    Some((i - ng as isize, j - ng as isize, d as usize))
}


/// Implementation of `Domain` for 2D cartesian meshes
impl <S: Clone + Default> Domain<CartesianBlock<S>> {
    pub fn new_rectangular(low: RealVec2, high: RealVec2,  n_cells: UIntVec2,
                           n_blocks: UIntVec2, var_names: Vec<String>,
                           num_comps: Vec<usize>, bcs: Vec<BCs<S>>,
                           num_ghost: Vec<usize> ) -> Domain<CartesianBlock<S>> {

        // the cells in the domain are divided evenly between each block
        // the floor of the division is taken in case the number of cells isn't divisible by
        // the number of blocks
        let total_blocks = n_blocks[0] * n_blocks[1];
        let n_cells_per_block_x = (n_cells[0])/(n_blocks[0]);
        let n_cells_per_block_y = (n_cells[1])/(n_blocks[1]);

        // the last block takes the left over cells
        let n_cells_final_x = n_cells[0] % n_blocks[0];
        let n_cells_final_y = n_cells[1] % n_blocks[1];

        // compute dx and dy to assist in finding the high and low corner of each block
        let dx = (high[0] - low[0])/(n_cells[0] as Real);
        let dy = (high[1] - low[1])/(n_cells[1] as Real);

        // construct each block
        let mut blocks: Vec<CartesianBlock<S>> = Vec::with_capacity(total_blocks);
        for i_block in 0..n_blocks[0] {
            for j_block in 0..n_blocks[1] {
                // the block's id is equivalent to its flattened index
                let id = get_flat_index(
                    i_block as isize, j_block as isize, 0, 
                    &n_blocks, 0
                ).unwrap();

                // compute the neighbouring blocks
                let north = get_flat_index(
                    i_block as isize, j_block as isize + 1, 0, &n_blocks, 0
                );
                let south = get_flat_index(
                    i_block as isize, j_block as isize - 1, 0, &n_blocks, 0
                );
                let east = get_flat_index(
                    i_block as isize + 1, j_block as isize, 0, &n_blocks, 0
                );
                let west = get_flat_index(
                    i_block as isize - 1, j_block as isize, 0, &n_blocks, 0
                );

                // compute the high corner and the number of cells for the block
                let mut cells = UIntVec2([n_cells_per_block_x, n_cells_per_block_y]);
                let mut high_block = RealVec2([
                    low[0] + (n_cells_per_block_x * (i_block + 1)) as Real * dx,
                    low[1] + (n_cells_per_block_y * (j_block + 1)) as Real * dy
                ]);
                if i_block == n_blocks[0] - 1{
                    cells[0] += n_cells_final_x;
                    high_block[0] += n_cells_final_x as Real * dx;
                }
                if j_block == n_blocks[1] - 1{
                    cells[1] += n_cells_final_y;
                    high_block[1] += n_cells_final_y as Real * dy;
                }

                // construct the block
                let block = CartesianBlock::<S>::new(
                    None,
                    RealVec2([
                        low[0] + (n_cells_per_block_x * i_block) as Real * dx,
                        low[1] + (n_cells_per_block_y * j_block) as Real * dy
                    ]),
                    high_block,
                    cells,
                    var_names.clone(),
                    num_comps.clone(),
                    num_ghost.clone(),
                    id,
                    // The boundary conditions need to be modified to account for internal boundaries
                    bcs.clone(), 
                    [west, south],
                    [east, north],
                );
                // push the block into the vector containing the blocks
                blocks.push(block);
            }
        }
        // construct and return the domain!
        Domain{blocks}
    }
}


#[cfg(test)]
mod tests{
    use super::*;

    #[test]
    fn rectangular_domain () {
        let _domain = Domain::<CartesianBlock<Real>>::new_rectangular(
            RealVec2([0.0, 0.0]),                               // low
            RealVec2([10.0, 10.0]),                             // high
            UIntVec2([5, 5]),                                   // n_cells
            UIntVec2([2,2]),                                    // n_blocks
            vec![String::from("var_1"), String::from("var_2")], // var_names
            vec![2, 1],                                         // n_comps
            vec![BCs::empty(), BCs::empty()],                   // bcd
            vec![1,1]                                           // num_ghost
        );
        //println!("domain = {:?}", _domain);
    }
}
