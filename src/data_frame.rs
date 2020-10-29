/// Module which defines and implements data types to store the data on a mesh
/// Currently, only CartesianMesh's are supported

use crate::mesh;

/// Available boundary conditions
#[allow(dead_code)]
pub enum BCType
{
    Prescribed(Vec<f64>),
    Dirichlet(f64),
    Neumann,
    // UseDefined (some function)
}

// ****************** Define traits all data frames must have *****************

/// BoundaryCondition trait for DataFrames
pub trait BoundaryCondition
{
    fn fill_bc(&mut self, bc_lo: &BCType, bc_hi: &BCType);
}

/// InitialCondition trait for DataFrames
pub trait InitialCondition 
{
    fn fill_ic (&mut self, ic: &dyn Fn(f64)->f64);
}




// *********************** Cartesian data frame *******************************

/// Structure to store data defined on a CartesianMesh
#[derive(Debug)]
pub struct CartesianDataFrame
{
    pub data: Vec<f64>,
    pub n_ghost: u32,
    pub underlying_mesh: mesh::CartesianMesh,
}

/// data structure to store data on CartesianMesh
impl CartesianDataFrame 
{
    /// generate new data structure from a given mesh, adding a given number of ghost nodes
    pub fn new_from(m: mesh::CartesianMesh, n_ghost: u32) -> CartesianDataFrame
    {
        let u = CartesianDataFrame
        {
            n_ghost,
            data: vec![0.0; m.n + 2*n_ghost as usize],
            underlying_mesh: m,
        };
        u
    }
}

// I don't know why this isn't accessible elsewhere if implementing CartesianDataFrame
/// Implement InitialCondition trait for CartesianDataFrame
impl CartesianDataFrame // for InitialCondition
{
    pub fn fill_ic (&mut self, ic: &dyn Fn(f64)->f64)
    {
        for (i,x) in self.underlying_mesh.node_pos.iter().enumerate()
        {
            self.data[i+self.n_ghost as usize] = ic(*x);
        }
    }
}

// I don't know why this isn't accessible elsewhere if implementing CartesianDataFrame
/// Implementation of boundary conditions for CartesianDataFrame
impl CartesianDataFrame // for Boundary Condition
{
    /// Fill the ghost nodes of the CartesianDataFrame based on BCType
    pub fn fill_bc (&mut self, bc_lo: &BCType, bc_hi: &BCType)
    {
        // low boundary condition
        match bc_lo
        {
            BCType::Prescribed(values) =>
            {
                let mut i = 0;
                for val in values.iter().rev()
                {
                    self.data[i] = *val;
                    i += 1;
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

/// test module for data frames. 
/// Currently tests implementation of both types of boundary condition
#[cfg(test)]
mod tests
{
    use crate::mesh;
    use crate::data_frame;
    use crate::BCType;

    fn initial_condition(x: f64) -> f64
    {
        x + 1.0
    }

    fn test_setup () -> data_frame::CartesianDataFrame
    {
        let u1 = mesh::CartesianMesh::new(0.0, 10.0, 5);
        let mut cdf = data_frame::CartesianDataFrame::new_from(u1, 2); 
        //note df now owns u1, and u1 is invalid

        cdf.fill_ic(&initial_condition);

        cdf
    }

    #[test]
    fn test_dirichlet ()
    {
        let mut cdf = test_setup();
        cdf.fill_bc(&BCType::Dirichlet(0.0), &BCType::Dirichlet(1.0));
        assert_eq!(cdf.data, vec![-6.0, -2.0, 2.0, 4.0, 6.0, 8.0, 10.0, -8.0, -26.0]);
    }

    #[test]
    fn test_neumann ()
    {
        let mut cdf = test_setup();
        cdf.fill_bc(&BCType::Neumann, &BCType::Neumann);
        assert_eq!(cdf.data, vec![4.0, 2.0, 2.0, 4.0, 6.0, 8.0, 10.0, 10.0, 8.0]);
    }
}