use super::mesh::*;
use super::dataframe::*;
use std::rc::Rc;


type DfHashMap<T> = std::collections::HashMap<String, CartesianDataFrame2D<T>>;


/// # Cartesian block
/// For parallel computations, the entire domain is split up into segments. A `CartesianBlock`
/// contains all the information for one of these segments, and where it sits in relation to the 
/// other segments
#[derive(Debug, Clone, Default)]
pub struct CartesianBlock<T> {
    /// An id to identify the block, and will be used for identifying neighbouring blocks
    pub id: usize,

    /// Optional name for the block, to make identifying it easier for people
    pub name: Option<String>,

    /// The underlying mesh for this block
    pub mesh: Rc<CartesianMesh2D>,

    /// The dataframes for the block, stored in a hashmap so they can be identified by a name
    pub dfs: DfHashMap<T>,

    /// The blocks stored on the low side of the block. None if there is no block
    pub lo_connectivity: [Option<usize>; 2],

    /// The blocks stored on the high side of the block. None if there is no block
    pub hi_connectivity: [Option<usize>; 2],
}



impl <T> core::ops::Index<&str> for CartesianBlock<T> {
    type Output = CartesianDataFrame2D<T>;

    fn index (&self, name: &str) -> &Self::Output {
        &self.dfs
             .get(name)
             .unwrap_or_else(|| 
                panic!("{} not in CartesianBlock with id: {}, name: {:?}", name, self.id, self.name)
            )
    }
}
