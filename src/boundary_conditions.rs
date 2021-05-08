use crate::Real;

// ****************** Define traits all data frames must have *****************
/// All data frames must have a `BoundaryCondition` trait, so that the ghost
/// nodes can be filled
pub trait BoundaryCondition {
    /// function which fills the boundary conditions
    fn fill_bc(&mut self, bcs: &BCs);

    /// checks if the cell at (i,j,k) contains a valid or a ghost cell. Returns true if valid,
    /// and returns false if ghost
    fn ijk_is_valid_cell(&self, i: isize, j: isize, k: isize) -> bool;

    /// check if the cell at p contains a valid or ghost cell. Returns same as ijk_is_valid_cell
    fn p_is_valid_cell(&self, p: usize) -> bool;
}

/// Available boundary conditions
#[allow(dead_code)]
#[derive(Debug, Clone, Copy)]
pub enum BcType {
    /// Linearly extrapolates ghost cells from the data in the valid region,
    /// so that the value on the boundary is the value supplied by the user
    Dirichlet(Real),

    /// Specifies the gradient at the boundary
    Neumann(Real),

    /// The ghost cells will become a reflection of the actual values
    Reflect,

    /// boundary is internal to the overall domain, and should be filled with values from the block
    /// with the given `id` 
    Internal(usize),
}

// /// specifies the boundary condition for each inner component of the data
// #[derive(Debug, Clone)]
// pub struct InnerCompBCs {
//     pub lo: Vec<BcType>,
//     pub hi: Vec<BcType>,
// }
// impl InnerCompBCs {
//     pub fn new(lo: Vec<BcType>, hi: Vec<BcType>)->InnerCompBCs{
//         InnerCompBCs{lo, hi,}
//     }
// }

/// Specifies the boundary condition for each component of the data
#[derive(Debug, Clone)]
pub struct ComponentBCs {
    /// Specifies the boundary condition for the lower edges
    pub lo: Vec<BcType>,

    /// Specifies the boundary condition for the upper edges
    pub hi: Vec<BcType>,
}
impl ComponentBCs {
    pub fn new(lo: Vec<BcType>, hi: Vec<BcType>) -> ComponentBCs{
        ComponentBCs {lo, hi,}
    }
}

/// Specifies the boundary conditions to apply for an entire data structure
#[derive(Debug, Clone)]
pub struct BCs {
    pub bcs: Vec<ComponentBCs>,
}
impl BCs {
    pub fn new(bcs: Vec<ComponentBCs>) -> BCs {
        BCs {
            bcs,
        }
    }

    pub fn empty() -> BCs{
        BCs{bcs: Vec::new()}
    }
}

