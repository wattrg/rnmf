// /// one dimensional mesh
// pub mod cartesian1d;

/// # Two dimensional cartesian mesh
/// A mesh where each cell is the same size, and has edges perpendicular to the *x* and *y* axis.
pub mod cartesian2d;

/// Main data structure describing the domain
pub struct Domain<M> {
    blocks: Vec<M>,
}
