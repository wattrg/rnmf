
// ****************** Define traits all data frames must have *****************
/// All data frames must have a `BoundaryCondition` trait, so that the ghost
/// nodes can be filled
pub trait BoundaryCondition {
    /// function which fills the boundary conditions
    fn fill_bc(&mut self, bc_lo: BCType, bc_hi: BCType);
}

/// Available boundary conditions
#[allow(dead_code)]
pub enum BCType {
    /// user supplies a vector which is placed into ghost cells.
    /// The first entry in the vector is put in the ghost cell closest to
    /// the valid computational domain
    Prescribed(Vec<f64>),

    /// Linearly extrapolates ghost cells from the data in the valid region,
    /// so that the value on the boundary is the value supplied by the user
    Dirichlet(f64),

    /// Evenly reflects data near boundary into ghost nodes
    Neumann,
}

