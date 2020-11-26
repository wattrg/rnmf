/// ## Module providing cartesian meshes
/// 
/// To create a new cartesian mesh, first bring the module into scope with 
/// `use crate::mesh::cartesian::*`. Then a new cartesian mesh is made in 2 steps. 
/// The first step is to create the nodes (this is referred to as the mesh): 
/// 
/// ``` rust
/// let mesh = CartesianMesh::new(lo: Vec<f64>, hi: Vec<f64> , n: Vec<int>, dim: usize);
/// ```
/// And then create one or more data structures to hold the data on top
/// of the mesh:
/// ```
/// let mut data1 = CartesianDataFrame::new_from(&mesh, n_comp: usize, n_ghost: usize);
/// let mut data2 = CartesianDataFrame::new_from(&mesh, n_comp: usize, n_ghost: usize);
/// ```
/// Following the creation of the data frame, the initial conditions can be filled by using the
/// `fill_ic` method, which accepts a function/closure:
/// ```rust
/// data1.fill_ic(&mut self, ic: fn(f64, f64, f64)->f64);
/// ```
/// Finally, the boundary conditions can be filled with the `fill_bc` method:
/// ```rust
/// data1.fill_bc(&mut self, bc_lo: BCType, bc_hi: BCType);
/// ```
/// Once this has been completed, the data is fully initialised, and is ready to be used.
pub mod cartesian;
