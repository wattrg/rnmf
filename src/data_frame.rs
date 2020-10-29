/// Module which defines and implements data types to store the data on a mesh
/// Currently, only CartesianMesh's are supported

use crate::mesh;

/// Available boundary conditions
#[allow(dead_code)]
pub enum BCType
{
    /// user supplies a vector which is placed into ghost cells.
    /// The first entry in the vector is put in the ghost cell closest to
    /// the valid computational domain
    Prescribed(Vec<f64>),

    /// Linearly extrapolates ghost cells from the data in the valid region,
    /// so that the value on the boundary is the value supplied by the user
    Dirichlet(f64),

    /// Evenly reflects data near boundary into ghost nodes
    Neumann,
    // UseDefined (some function)
}

// ****************** Define traits all data frames must have *****************

/// All data frames must have a `BoundaryCondition` trait, so that the ghost
/// nodes can be filled
pub trait BoundaryCondition
{
    /// function which fills the boundary conditions
    fn fill_bc(&mut self, bc_lo: &BCType, bc_hi: &BCType);
}



// *********************** Cartesian data frame *******************************

/// Structure to store data defined on a `CartesianMesh`
#[derive(Debug)]
pub struct CartesianDataFrame <'a>
{
    /// The data is stored here
    pub data: Vec<f64>,

    /// The number of ghost nodes, added to each side of the underlying `CartesianMesh`
    pub n_ghost: u32,

    /// Reference to the underlying `CartesianMesh`
    pub underlying_mesh: &'a mesh::CartesianMesh,
}

/// data structure to store data on CartesianMesh
impl CartesianDataFrame <'_>
{
    /// generate new `CartesianDataFrame` from a given mesh, adding a given 
    /// number of ghost nodes
    pub fn new_from(m: & mesh::CartesianMesh, n_ghost: u32) -> CartesianDataFrame
    {
        CartesianDataFrame
        {
            n_ghost,
            data: vec![0.0; m.n + 2*n_ghost as usize],
            underlying_mesh: m,
        }
    }


    /// Fill `CartesianDataFrame` from a initial condition function
    pub fn fill_ic (&mut self, ic: &dyn Fn(f64)->f64)
    {
        for (i,x) in self.underlying_mesh.node_pos.iter().enumerate()
        {
            self.data[i+self.n_ghost as usize] = ic(*x);
        }
    }
}

/// Implementation of `BoundaryCondition` for `CartesianDataFrame`
impl BoundaryCondition for CartesianDataFrame <'_> // for Boundary Condition
{
    /// Fill the ghost nodes of the CartesianDataFrame based on BCType
    fn fill_bc (&mut self, bc_lo: &BCType, bc_hi: &BCType)
    {
        // low boundary condition
        match bc_lo
        {
            BCType::Prescribed(values) =>
            {
                //let mut i = 0;
                for (i,val) in values.iter().rev().enumerate()
                {
                    self.data[i] = *val;
                    //i += 1;
                }

            }
            BCType::Neumann =>
            {
                for i in 0usize..self.n_ghost as usize
                {
                    self.data[self.n_ghost as usize - i - 1] = self.data[self.n_ghost as usize + i];
                }
            }
            BCType::Dirichlet(value) =>
            {
                let m: f64 = self.data[self.n_ghost as usize] - value;
                for i in 0usize..self.n_ghost as usize
                {
                    self.data[self.n_ghost as usize - i - 1] = self.data[self.n_ghost as usize]-2.0*(i as f64+1.0)*m
                }
            }
        }

        // high boundary condition
        let n: usize = self.data.len()-1;
        match bc_hi
        {
            BCType::Prescribed(values) =>
            {
                for (i, val) in values.iter().enumerate()
                {
                    self.data[i+self.n_ghost as usize + self.underlying_mesh.n] = *val;
                }
            }
            BCType::Neumann =>
            {
                for i in 1usize..self.n_ghost as usize + 1
                {
                    self.data[n-(self.n_ghost as usize)+i] = self.data[n-self.n_ghost as usize-i+1];
                }
            }
            BCType::Dirichlet(value) =>
            {
                let m = value - self.data[n-self.n_ghost as usize];
                for i in 1usize..self.n_ghost as usize + 1
                {
                    self.data[n-(self.n_ghost as usize) + i] = self.data[n-self.n_ghost as usize] + 2.0*m*i as f64;
                }
            }
        }
    }
}








// test module for data frames. 
// Currently tests implementation of boundary conditions
#[cfg(test)]
mod tests
{
    use crate::mesh;
    use crate::data_frame;
    use crate::{BCType, BoundaryCondition};

    fn initial_condition(x: f64) -> f64
    {
        x + 1.0
    }

    #[test]
    fn test_dirichlet_bc ()
    {
        let u1 = mesh::CartesianMesh::new(0.0, 10.0, 5);
        let mut cdf = data_frame::CartesianDataFrame::new_from(&u1, 2);
        cdf.fill_ic(&initial_condition);
        cdf.fill_bc(&BCType::Dirichlet(0.0), &BCType::Dirichlet(1.0));
        assert_eq!(cdf.data, vec![-6.0, -2.0, 2.0, 4.0, 6.0, 8.0, 10.0, -8.0, -26.0]);
    }

    #[test]
    fn test_neumann_bc ()
    {
        let u1 = mesh::CartesianMesh::new(0.0, 10.0, 5);
        let mut cdf = data_frame::CartesianDataFrame::new_from(&u1, 2);
        cdf.fill_ic(&initial_condition);
        cdf.fill_bc(&BCType::Neumann, &BCType::Neumann);
        assert_eq!(cdf.data, vec![4.0, 2.0, 2.0, 4.0, 6.0, 8.0, 10.0, 10.0, 8.0]);
    }

    #[test]
    fn test_prescribed_bc()
    {
        let u1 = mesh::CartesianMesh::new(0.0, 10.0, 5);
        let mut cdf = data_frame::CartesianDataFrame::new_from(&u1, 2);
        cdf.fill_ic(&initial_condition);
        cdf.fill_bc(&BCType::Prescribed(vec![-1.0, -2.0]), &BCType::Prescribed(vec![15.0, 16.0]));
        assert_eq!(cdf.data, vec![-2.0, -1.0, 2.0, 4.0, 6.0, 8.0, 10.0, 15.0, 16.0]);
    }
}