// ****************** Define traits all data frames must have *****************
/// All data frames must have a `BoundaryCondition` trait, so that the ghost
/// nodes can be filled
pub trait BoundaryCondition <S> {
    /// function which fills the boundary conditions
    fn fill_bc(&mut self);
}

/// Available boundary conditions
#[derive(Debug, Clone)]
pub enum BcType <S> {
    /// Linearly extrapolates ghost cells from the data in the valid region,
    /// so that the value on the boundary is the value supplied by the user
    Dirichlet(S),

    /// Specifies the gradient at the boundary
    Neumann(S),

    /// The ghost cells will become a reflection of the actual values
    Reflect,

    /// The values in the ghost cells are prescribed by a vector of values
    Prescribed(Vec<S>),
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
pub struct ComponentBCs<S> {
    /// Specifies the boundary condition for the lower edges
    pub lo: Vec<BcType<S>>,

    /// Specifies the boundary condition for the upper edges
    pub hi: Vec<BcType<S>>,
}

impl <S> ComponentBCs<S> {
    pub fn new(lo: Vec<BcType<S>>, hi: Vec<BcType<S>>) -> ComponentBCs<S>{
        ComponentBCs {lo, hi,}
    }
}

/// Specifies the boundary conditions to apply for an entire data structure
#[derive(Debug, Clone)]
pub struct BCs<S> {
    pub bcs: Vec<ComponentBCs<S>>,
}

impl <S> BCs<S> {
    pub fn new(bcs: Vec<ComponentBCs<S>>) -> BCs<S> {
        BCs {
            bcs,
        }
    }

    pub fn empty() -> BCs<S>{
        BCs{bcs: Vec::new()}
    }
}

