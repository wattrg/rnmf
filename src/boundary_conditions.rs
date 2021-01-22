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
pub enum BCType {
    /// user supplies a vector which is placed into ghost cells.
    /// The first entry in the vector is put in the ghost cell closest to
    /// the valid computational domain
    Prescribed(Vec<Real>),

    /// Linearly extrapolates ghost cells from the data in the valid region,
    /// so that the value on the boundary is the value supplied by the user
    Dirichlet(Real),

    /// Specifies the gradient at the boundary
    Neumann(Real),
}

/// Specifies the boundary condition for each component of the data
pub struct ComponentBCs {
    /// Specifies the boundary condition for the lower edges
    pub lo: Vec<BCType>,

    /// Specifies the boundary condition for the upper edges
    pub hi: Vec<BCType>,
}


impl ComponentBCs {
    pub fn new(lo: Vec<BCType>, hi: Vec<BCType>) -> ComponentBCs{
        ComponentBCs {
            lo,
            hi,
        }
    }
}

/// Specifies the boundary conditions to apply for an entire data structure
pub struct BCs {
    pub bcs: Vec<ComponentBCs>,
}
impl BCs {
    pub fn new(bcs: Vec<ComponentBCs>) -> BCs {
        BCs {
            bcs,
        }
    }
}

