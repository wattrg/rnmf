/// Module which defines and implements data types to store the data on a mesh
/// Currently, only CartesianMesh's are supported
use crate::mesh;
use crate::boundary_conditions::{BCType, BoundaryCondition};





// *********************** Cartesian data frame *******************************
/// Structure to store data defined on a `CartesianMesh`
#[derive(Debug)]
pub struct CartesianDataFrame<'a>{
    /// The data is stored here
    pub data: Vec<f64>,

    /// The number of ghost nodes, added to each side of the underlying `CartesianMesh`
    pub n_ghost: usize,

    n_grown: Vec<usize>,

    /// Reference to the underlying `CartesianMesh`
    pub underlying_mesh: &'a mesh::CartesianMesh,

    pub n_comp: usize,

    n_nodes: usize,
}

/// data structure to store data on CartesianMesh
impl CartesianDataFrame<'_> {
    /// generate new `CartesianDataFrame` from a given mesh, adding a given
    /// number of ghost nodes
    pub fn new_from(m: & mesh::CartesianMesh, 
                    n_comp: usize, 
                    n_ghost: usize) -> CartesianDataFrame
    {
        // this could be generated starting from the number of nodes in the 
        // underlying mesh, but I think this is just as fast in that it requires
        // the same number of operations
        let mut n_nodes = m.n[0] + 2*n_ghost;
        if m.dim >= 2
        {
            n_nodes *= m.n[1] + 2*n_ghost;

            if m.dim == 3
            {
                n_nodes *= m.n[2] + 2*n_ghost;
            }
        }



        CartesianDataFrame
        {
            n_ghost,
            data:  vec![0.0; n_nodes],
            underlying_mesh: m,
            n_comp,
            n_nodes,
            n_grown: m.n.clone().into_iter().map(|n| n + 2*n_ghost).collect(),
        }
    }

    /// Fill `CartesianDataFrame` from a initial condition function
    pub fn fill_ic (&mut self, ic: fn(f64, f64, f64)->f64)
    {
        for (i,x) in self.underlying_mesh.node_pos.iter().enumerate()
        {
            self.data[i+self.n_ghost] = ic(*x, 0.0, 0.0);
        }
    }


    // unused for the moment, but will form the basis for indexing the data frames
    /// Retrieves the element at (i,j,k,n). The valid cells are index from zero
    /// and ghost cells at the lower side of the domain are indexed with negative
    /// numbers.
    #[allow(dead_code)] 
    fn index_mut(&mut self, i: usize, j: usize, k: usize, n: usize) -> f64
    {
        self.data[n + 
                  self.n_comp*(i+self.n_ghost) +
                  self.underlying_mesh.n[1]*self.n_comp*(j+self.n_ghost)+
                  self.n_comp*self.underlying_mesh.n[1]*self.underlying_mesh.n[2]*(k+self.n_ghost)]
    }


}




/// Implementation of `BoundaryCondition` for `CartesianDataFrame`
impl BoundaryCondition for CartesianDataFrame<'_> {
    /// Fill the ghost nodes of the CartesianDataFrame based on BCType
    fn fill_bc(&mut self, bc_lo: BCType, bc_hi: BCType) {
        // low boundary condition
        match bc_lo {
            BCType::Prescribed(values) => {
                for (i, &val) in values.iter().rev().enumerate() {
                    self.data[i] = val;
                }
            }
            BCType::Neumann => {
                for i in 0usize..self.n_ghost as usize {
                    self.data[self.n_ghost - i - 1] = self.data[self.n_ghost + i];
                }
            }
            BCType::Dirichlet(value) => {
                let m: f64 = self.data[self.n_ghost as usize] - value;
                for i in 0usize..self.n_ghost as usize {
                    self.data[self.n_ghost - i - 1] =
                        self.data[self.n_ghost] - 2.0 * (i as f64 + 1.0) * m
                }
            }
        }

        // high boundary condition
        let n: usize = self.data.len() - 1;
        match bc_hi {
            BCType::Prescribed(values) => {
                for (i, val) in values.iter().enumerate() {
                    self.data[i + self.n_ghost as usize + self.underlying_mesh.n[0]] = *val;
                }
            }
            BCType::Neumann => {
                for i in 1usize..self.n_ghost as usize + 1 {
                    self.data[n - (self.n_ghost as usize) + i] =
                        self.data[n - self.n_ghost as usize - i + 1];
                }
            }
            BCType::Dirichlet(value) => {
                let m = value - self.data[n - self.n_ghost as usize];
                for i in 1usize..self.n_ghost as usize + 1 {
                    self.data[n - (self.n_ghost as usize) + i] =
                        self.data[n - self.n_ghost as usize] + 2.0 * m * i as f64;
                }
            }
        }
    }
}

// test module for data frames.
// Currently tests implementation of boundary conditions
#[cfg(test)]
mod tests {
    use crate::data_frame;
    use crate::mesh;
    use crate::{BCType, BoundaryCondition};

    #[test]
    fn test_dirichlet_bc() {
        let u1 = mesh::CartesianMesh::new(vec![0.0], vec![10.0], vec![5], 1);
        let mut cdf = data_frame::CartesianDataFrame::new_from(&u1, 1, 2);
        cdf.fill_ic(|x,_,_| x + 1.0);
        cdf.fill_bc(BCType::Dirichlet(0.0), BCType::Dirichlet(1.0));
        assert_eq!(
            cdf.data,
            vec![-6.0, -2.0, 2.0, 4.0, 6.0, 8.0, 10.0, -8.0, -26.0]
        );
    }

    #[test]
    fn test_neumann_bc() {
        let u1 = mesh::CartesianMesh::new(vec![0.0], vec![10.0], vec![5], 1);
        let mut cdf = data_frame::CartesianDataFrame::new_from(&u1, 1, 2);
        cdf.fill_ic(|x, _,_| x + 1.0);
        cdf.fill_bc(BCType::Neumann, BCType::Neumann);
        assert_eq!(
            cdf.data,
            vec![4.0, 2.0, 2.0, 4.0, 6.0, 8.0, 10.0, 10.0, 8.0]
        );
    }

    #[test]
    fn test_prescribed_bc() {
        let u1 = mesh::CartesianMesh::new(vec![0.0], vec![10.0], vec![5], 1);
        let mut cdf = data_frame::CartesianDataFrame::new_from(&u1, 1, 2);
        cdf.fill_ic(|x, _,_| x + 1.0);
        cdf.fill_bc(
            BCType::Prescribed(vec![-1.0, -2.0]),
            BCType::Prescribed(vec![15.0, 16.0]),
        );
        assert_eq!(
            cdf.data,
            vec![-2.0, -1.0, 2.0, 4.0, 6.0, 8.0, 10.0, 15.0, 16.0]
        );
    }
}
