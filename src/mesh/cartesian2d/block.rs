use super::mesh::*;
use super::dataframe::*;
use std::rc::Rc;
use super::super::Block;
use crate::{RealVec2, UIntVec2};
use crate::boundary_conditions::BCs;


/// # Storage for the `CartesianDataFrame2D` associated with a block
/// In case multiple data frames are needed for a given `Block` (e.g. if you have multiple 
/// variables which don't logically fit together, such as two fluid species or fluids and 
/// electromagnetic fields), the dataframes are stored behind a hashmap, allowing them to be 
/// accessed by a name easily recognisable by humans.
type DfHashMap<T> = std::collections::HashMap<String, CartesianDataFrame2D<T>>;

type Connectivity = [Option<usize>; 2];

/// # Cartesian block
/// For parallel computations, the entire domain is split up into segments. A `CartesianBlock`
/// contains all the information for one of these segments, and where it sits in relation to the 
/// other segments
#[derive(Debug)]
pub struct CartesianBlock<T: Clone + Default> {
    /// An id to identify the block, and will be used for identifying neighbouring blocks
    pub id: usize,

    /// Optional name for the block, to make identifying it easier for people
    pub name: Option<String>,

    /// The underlying mesh for this block
    mesh: Rc<CartesianMesh2D>,

    /// The dataframes for the block, stored in a hashmap so they can be identified by a name
    dfs: DfHashMap<T>,

    /// The blocks stored on the low side of the block. None if there is no block
    low_connect: Connectivity,

    /// The blocks stored on the high side of the block. None if there is no block
    high_connect: Connectivity,
}

impl <T: Clone + Default> Block for CartesianBlock<T>{}

impl <T: Clone+Default> core::ops::Index<String> for CartesianBlock<T> {
    type Output = CartesianDataFrame2D<T>;

    fn index (&self, name: String) -> &Self::Output {
        &self.dfs
             .get(&name)
             .unwrap_or_else(|| 
                panic!("{} not in CartesianBlock with id: {}, name: {:?}", name, self.id, self.name)
            )
    }
}

impl <T: Clone + Default> CartesianBlock<T> {
    pub fn new(name: Option<String>, low: RealVec2, high: RealVec2,
               n_cells: UIntVec2, data_names: Vec<String>, num_comps: Vec<usize>,
               num_ghost: Vec<usize>, id: usize, bcs: Vec<BCs>,
               low_connect: Connectivity, high_connect: Connectivity) -> CartesianBlock<T>{

        // create the mesh
        let mesh = CartesianMesh2D::new(low, high, n_cells);

        // create the dataframes
        let mut dfs = DfHashMap::<T>::new();
        for (i, name) in data_names.iter().enumerate() {
            let df = CartesianDataFrame2D::<T>::new_from(
                &mesh, bcs[i].clone(), num_comps[i], num_ghost[i]
            );
            dfs.insert(name.to_string(), df);
        }

        // create the block
        CartesianBlock{
            name,
            id,
            dfs,
            mesh,
            low_connect,
            high_connect,
        }
    }
}



#[cfg(test)]
mod tests {
    use super::*;
    use crate::Real;

    #[test]
    fn new_cartesian_block() {
        let cb = CartesianBlock::<Real>::new(
            Some(String::from("test")),
            RealVec2([5.0, 5.0]),
            RealVec2([10.0, 10.0]),
            UIntVec2([5, 5]),
            vec![String::from("var_1"), String::from("var_2")],
            vec![1, 1],
            vec![1, 1],
            0,
            vec![BCs::empty(), BCs::empty()],
            [None, None],
            [None, None],
        );

        //println!("{:?}", cb);
    }
}

