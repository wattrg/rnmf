
use std::rc::Rc;


/// Structure containing data to define a CartesianMesh
#[derive(Debug)]
pub struct CartesianMesh {
    /// The position of the nodes in the mesh
    pub node_pos: Vec<f64>,

    /// The lower corner of the mesh
    pub lo: Vec<f64>,

    /// The upper corner of the mesh
    pub hi: Vec<f64>,

    /// The number of cells in nodes in each direction
    pub n : Vec<usize>,

    /// The distance between each node in each
    pub dx: Vec<f64>,

    /// The number of spatial dimensions of the mesh
    pub dim: usize,
    
    /// The number of nodes in the mesh
    n_nodes: usize
}

trait Indexing{
    /// Returns the non-flat index of the item stored at a particular flattened location
    fn get_ijk(&self, p: usize ) -> Option<(usize, usize, usize, usize)>;

    /// Returns the value stored at a particular non-flat index
    fn index(&self, i: usize, j: usize, k: usize, n: usize) -> f64;
}


impl CartesianMesh {
    /// Generate new CartesianMesh from the lo and high corner, 
    /// and the number of nodes
    pub fn new(lo: Vec<f64>, hi: Vec<f64>, n:Vec<usize>, dim: usize) -> Rc<CartesianMesh>
    {
        let mut dx = vec![0.0; dim as usize];
        for i_dim in 0..dim as usize{
            dx[i_dim] = (hi[i_dim] - lo[i_dim])/(n[i_dim] as f64);
        }

        // calculate the capacity of the vector required to store all the node positions
        let mut n_nodes = n[0];
        if dim >= 2{
            n_nodes *= n[1];

            if dim == 3 {
                n_nodes *= n[2];
            }
        }
        n_nodes *= dim;

        // allocate memory to the node positions vector
        let node_pos = Vec::with_capacity(n_nodes);

        // allocate memory to the mesh
        let mut cm = CartesianMesh
        {
            lo,
            hi,
            n,
            dx,
            node_pos,
            dim,
            n_nodes,
        };

        // calculate the positions of the nodes
        cm.node_pos = (0..n_nodes).map(
            |p: usize| -> f64 {
                let (i,j,k,n) = cm.get_ijk(p).unwrap();
                match dim {
                    1 => {
                        cm.lo[n] + ((i as f64) + 0.5) * cm.dx[n]
                    }
                    2 => {
                        match n {
                            0 => {cm.lo[n] + ((i as f64) + 0.5) * cm.dx[n]}
                            1 => {cm.lo[n] + ((j as f64) + 0.5) * cm.dx[n]}
                            _ => {panic!("n cannot be {} in {}D", n, dim)}
                        }
                    }
                    3 => {
                        match n{
                            0 => {cm.lo[n] + ((i as f64) + 0.5) * cm.dx[n]}
                            1 => {cm.lo[n] + ((j as f64) + 0.5) * cm.dx[n]}
                            2 => {cm.lo[n] + ((k as f64) + 0.5) * cm.dx[n]}
                            _ => {panic!("n cannot be {} in {}D", n, dim)}
                        }
                    }
                    _ => {
                        panic!("{}D not supported!", dim);
                    }
                }
                
            }).collect();
        Rc::new(cm)
    }
}



impl Indexing for CartesianMesh
{
    /// Retrieves the element at (i,j,k,n)
    #[allow(dead_code)]
    fn index(&self, i: usize, j: usize, k: usize, n: usize) -> f64
    {
        let mut p = n + self.dim * i;
        if self.dim >= 2 {
            p += self.dim*self.n[1]*j;
        
            if self.dim == 3 {
                p += self.dim*self.n[2]*self.n[1]*k;
            }
        }
        
        let np = self.node_pos.get(p);
        match np{
            Some(np) => {*np},
            None => panic!("Index ({}, {}, {}, {}) out of range of Cartesian mesh with size {:?}", 
                                i,j,k,n, self.n),
        }

    }

    /// Returns the un-flattened index
    fn get_ijk(& self, p: usize) -> Option<(usize, usize, usize, usize)>{
        if p >= self.n_nodes {
            None
        }
        else{
            // the idea with the closures is that I'll implement lazy evaluation for each 
            // closure following https://doc.rust-lang.org/book/ch13-01-closures.html, so that we
            // aren't performing any unnecessary or repeated calculations
            let d = p % self.dim;
            
            let i = || {
                ((p-d)/self.dim) % self.n[0]
            };
            let j = || {
                ((p - d - i()*self.dim)/(self.dim * self.n[0])) % self.n[1]
            };
            let k = || {
                (p- j()*self.n[0]*self.dim - i()*self.dim - d)/(self.dim * self.n[0] * self.n[1])
            };
            

            match self.dim{
                1 => {
                    Some((i(), 0, 0, d))
                }
                2 => {
                    Some((i(), j(), 0, d))
                }
                3 => {
                    Some((i(), j(), k(), d))
                }
                _ => {
                    panic!("{}D not supported!", self.dim);
                }
            }
        }
    }
}


/// Immutable Iterator struct for CartesianMesh
pub struct CartesianMeshIter <'a> {
    mesh: &'a CartesianMesh,
    current_indx: usize,
}

impl <'a> Iterator for CartesianMeshIter<'a>{
    // return a vector of references instead of slice because the length of the returned object
    // depends on the number of dimensions, and thus is not known at compile time, so slices
    // won't work
    type Item = Vec<&'a f64>;
    
    fn next(& mut self) -> Option<Self::Item>{
        // If the next position doesn't exist, return None
        if self.current_indx > self.mesh.node_pos.len() -self.mesh.dim {
            None
        }
        // If the next position exists, return a vector of size self.dim containing references
        // to the position
        else {
            let mut next_pos = Vec::with_capacity(self.mesh.dim);
            for i_dim in 0..(self.mesh.dim) {
                next_pos.push(& self.mesh.node_pos[self.current_indx + i_dim])
            }
            self.current_indx += self.mesh.dim;
            Some(next_pos)
        }
    }
}

impl <'a> IntoIterator for &'a CartesianMesh {
    type Item = Vec<&'a f64>;
    type IntoIter = CartesianMeshIter<'a>;

    fn into_iter(self) -> CartesianMeshIter<'a> {
        CartesianMeshIter{
            current_indx: 0,
            mesh: & self,
        }
    }
}


use crate::boundary_conditions::{BCType, BoundaryCondition};

// ********************************* Cartesian data frame ****************************************
/// Structure to store data defined on a `CartesianMesh`
#[derive(Debug)]
pub struct CartesianDataFrame{
    /// The data is stored here
    pub data: Vec<f64>,

    /// The number of ghost nodes, added to each side of the underlying `CartesianMesh`
    pub n_ghost: usize,

    n_grown: Vec<usize>,

    /// Reference to the underlying `CartesianMesh`
    pub underlying_mesh: Rc<CartesianMesh>,

    pub n_comp: usize,

    n_nodes: usize,
}

/// data structure to store data on CartesianMesh
impl CartesianDataFrame {
    /// generate new `CartesianDataFrame` from a given mesh, adding a given
    /// number of ghost nodes
    pub fn new_from(m: & Rc<CartesianMesh>, 
                    n_comp: usize, n_ghost: usize) -> CartesianDataFrame
    {
        // this could be generated starting from the number of nodes in the 
        // underlying mesh, but I think this is just as fast in that it requires
        // the same number of operations
        let mut n_nodes = m.n[0] + 2*n_ghost;
        if m.dim >= 2 {
            n_nodes *= m.n[1] + 2*n_ghost;

            if m.dim == 3 {
                n_nodes *= m.n[2] + 2*n_ghost;
            }
        }
        n_nodes *= m.dim;


        CartesianDataFrame
        {
            n_ghost,
            data:  vec![0.0; n_nodes],
            n_grown: m.n.clone().into_iter().map(|n| n + 2*n_ghost).collect(),
            underlying_mesh: Rc::clone(m),
            n_comp,
            n_nodes,
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
}


impl Indexing for CartesianDataFrame{
    #[allow(dead_code)] 
    /// Retrieves the element at (i,j,k,n). The valid cells are index from zero
    /// and ghost cells at the lower side of the domain are indexed with negative
    /// numbers.
    fn index(&self, i: usize, j: usize, k: usize, n: usize) -> f64
    {
        self.data[n + 
                  self.n_comp*(i+self.n_ghost) +
                  self.underlying_mesh.n[1]*self.n_comp*(j+self.n_ghost)+
                  self.n_comp*self.underlying_mesh.n[1]*self.underlying_mesh.n[2]*(k+self.n_ghost)]
    }

    fn get_ijk(&self, p: usize) -> Option<(usize, usize, usize, usize)> {
        if p >= self.n_nodes {
            None
        }
        else{
            // the idea with the closures is that I'll implement lazy evaluation for each 
            // closure following https://doc.rust-lang.org/book/ch13-01-closures.html, so that we
            // aren't performing any unnecessary or repeated calculations
            let n0 = self.underlying_mesh.n[0];
            let n1 = self.underlying_mesh.n[1];
            let n = p % self.n_comp;
            
            let i = || {
                (p-n/self.n_comp) % n0
            };
            let j = || {
                ((p - n - i()*self.n_comp)/(self.n_comp * n0)) % n1
            };
            let k = || {
                (p- j()*n0*self.n_comp - i()*self.n_comp - n)/(self.n_comp * n0 * n1)
            };
            

            match self.underlying_mesh.dim{
                1 => {
                    Some((i() - self.n_ghost, 0, 0, n))
                }
                2 => {
                    Some((i() - self.n_ghost, j() - self.n_ghost, 0, n))
                }
                3 => {
                    Some((i() - self.n_ghost, j() - self.n_ghost, k() - self.n_ghost, n))
                }
                _ => {
                    panic!("{}D not supported!", self.underlying_mesh.dim);
                }
            }
        }
    }
}




impl <'a> BoundaryCondition for CartesianDataFrame{
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

/// mutable iterator for `CartesianDataFrame`.
pub struct CartesianDataFrameIter<'a> {
    pub df: &'a mut CartesianDataFrame,
    pub current_indx: usize,
}

impl <'a> Iterator for CartesianDataFrameIter<'a> {
    type Item = &'a mut f64;

    // this is a safe function containing unsafe code, since rust doesn't allow multiple mutable
    // references to independent elements of a vector. But in this case, because we are ever only
    // accessing each item once, it is actually safe.
    fn next(&mut self) -> Option<Self::Item> {
        // check if the next item exists
        if self.current_indx <= self.df.n_nodes-self.df.n_comp{
            // create raw pointer to the data
            let ptr = self.df.data.as_mut_ptr();
            // create variable to store the next item
            let next_data: Option<Self::Item>;
            unsafe{
                // access next item
                next_data = ptr.add(self.current_indx).as_mut();
            }
            // increment current index
            self.current_indx += self.df.n_comp;

            // return the next item
            next_data
        }
        // if the next item doesn't exist, return Option::None
        else {
            None
        }
    }
}

impl <'a> IntoIterator for &'a mut CartesianDataFrame{
    type Item = &'a mut f64;
    type IntoIter = CartesianDataFrameIter<'a>;

    fn into_iter(self) -> Self::IntoIter{
        Self::IntoIter {
            current_indx: 0,
            df: self,
        }
    }
}




#[cfg(test)]
mod tests{
    use super::*;

    #[test]
    // currently only tests the iterator on 1D data
    fn data_frame_iterator() {
        let m1 = CartesianMesh::new(vec![0.0], vec![6.0], vec![3], 1);
        let mut df = CartesianDataFrame::new_from(& m1, 1, 0);
        df.fill_ic(|x,_,_| x + 1.0);

        let mut df_iter = df.into_iter();
        assert_eq!(df_iter.next(), Some(&mut 2.0));
        assert_eq!(df_iter.next(), Some(&mut 4.0));
        assert_eq!(df_iter.next(), Some(&mut 6.0));
        assert_eq!(df_iter.next(), Option::None);
        
        // test to make sure mutating the data works
        for data in &mut df{
            *data += 1.0;
        }

        assert_eq!(df.data, vec![3.0, 5.0, 7.0]);
    }

    #[test]
    fn mesh_iterator () {
        let m = CartesianMesh::new(vec![0.0], vec![10.0], vec![5], 1);
        let mut m_iter = m.into_iter();
        assert_eq!(m_iter.next().unwrap(), vec![&mut 1.0]);
        assert_eq!(m_iter.next().unwrap(), vec![&mut 3.0]);
        assert_eq!(m_iter.next().unwrap(), vec![&mut 5.0]);
        assert_eq!(m_iter.next().unwrap(), vec![&mut 7.0]);
        assert_eq!(m_iter.next().unwrap(), vec![&mut 9.0]);
        assert_eq!(m_iter.next(), Option::None);

        let m2 = CartesianMesh::new(vec![0.0, 0.0], vec![6.0, 6.0], vec![3, 3], 2);
        let mut m2_iter = m2.into_iter();
        assert_eq!(m2_iter.next().unwrap(), vec![&mut 1.0, &mut 1.0]);
        assert_eq!(m2_iter.next().unwrap(), vec![&mut 3.0, &mut 1.0]);
        assert_eq!(m2_iter.next().unwrap(), vec![&mut 5.0, &mut 1.0]);
        assert_eq!(m2_iter.next().unwrap(), vec![&mut 1.0, &mut 3.0]);
        assert_eq!(m2_iter.next().unwrap(), vec![&mut 3.0, &mut 3.0]);
        assert_eq!(m2_iter.next().unwrap(), vec![&mut 5.0, &mut 3.0]);
        assert_eq!(m2_iter.next().unwrap(), vec![&mut 1.0, &mut 5.0]);
        assert_eq!(m2_iter.next().unwrap(), vec![&mut 3.0, &mut 5.0]);
        assert_eq!(m2_iter.next().unwrap(), vec![&mut 5.0, &mut 5.0]);
        assert_eq!(m_iter.next(), Option::None);
    }

    #[test]
    fn indexing () {
        // 1D cartesian mesh
        let m1 = CartesianMesh::new(vec![0.0], vec![10.0], vec![5], 1);
        assert_eq!(m1.get_ijk(3).unwrap(), (3,0,0,0));
        assert_eq!(m1.get_ijk(5), Option::None);
        assert_eq!(m1.index(3, 0, 0, 0), 7.0);
        let mut cdf1 = CartesianDataFrame::new_from(&m1, 1, 2);
        cdf1.fill_ic(|x,_,_|x + 1.0);

        let m2 = CartesianMesh::new(vec![0.0,0.0], vec![10.0,8.0], vec![5,4], 2);
        assert_eq!(m2.get_ijk(3).unwrap(), (1,0,0,1));
        assert_eq!(m2.get_ijk(15).unwrap(), (2,1,0,1));
        assert_eq!(m2.get_ijk(40), Option::None);
        assert_eq!(m2.index(2,1,0,1), 3.0);

        let m3 = CartesianMesh::new(vec![0.0,0.0,0.0], vec![10.0,8.0,10.0], vec![5,4,5], 3);
        assert_eq!(m3.get_ijk(3).unwrap(), (1,0,0,0));
        assert_eq!(m3.get_ijk(22).unwrap(), (2,1,0,1));
        assert_eq!(m3.get_ijk(60).unwrap(), (0,0,1,0));
        assert_eq!(m3.get_ijk(300), Option::None);
        assert_eq!(m3.index(2,1,0,1), 3.0);
    }

    #[test]
    fn node_pos () {
        let m1 = CartesianMesh::new(vec![0.0], vec![10.0], vec![5], 1);
        assert_eq!(m1.node_pos, vec![1.0, 3.0, 5.0, 7.0, 9.0]);

        let m2 = CartesianMesh::new(vec![0.0, 0.0], vec![10.0, 10.0], vec![5,5], 2);
        assert_eq!(
            m2.node_pos,
            vec![
                1.0, 1.0, 3.0, 1.0, 5.0, 1.0, 7.0, 1.0, 9.0, 1.0,
                1.0, 3.0, 3.0, 3.0, 5.0, 3.0, 7.0, 3.0, 9.0, 3.0,
                1.0, 5.0, 3.0, 5.0, 5.0, 5.0, 7.0, 5.0, 9.0, 5.0,
                1.0, 7.0, 3.0, 7.0, 5.0, 7.0, 7.0, 7.0, 9.0, 7.0,
                1.0, 9.0, 3.0, 9.0, 5.0, 9.0, 7.0, 9.0, 9.0, 9.0
            ]
        );

        let m3 = CartesianMesh::new(vec![0.0,0.0,0.0], vec![6.0,6.0,6.0], vec![3,3,3], 3);
        assert_eq!(
            m3.node_pos,
            vec![
                1.0, 1.0, 1.0, 3.0, 1.0, 1.0, 5.0, 1.0, 1.0,
                1.0, 3.0, 1.0, 3.0, 3.0, 1.0, 5.0, 3.0, 1.0,
                1.0, 5.0, 1.0, 3.0, 5.0, 1.0, 5.0, 5.0, 1.0,
                //
                1.0, 1.0, 3.0, 3.0, 1.0, 3.0, 5.0, 1.0, 3.0,
                1.0, 3.0, 3.0, 3.0, 3.0, 3.0, 5.0, 3.0, 3.0,
                1.0, 5.0, 3.0, 3.0, 5.0, 3.0, 5.0, 5.0, 3.0,
                //
                1.0, 1.0, 5.0, 3.0, 1.0, 5.0, 5.0, 1.0, 5.0,
                1.0, 3.0, 5.0, 3.0, 3.0, 5.0, 5.0, 3.0, 5.0,
                1.0, 5.0, 5.0, 3.0, 5.0, 5.0, 5.0, 5.0, 5.0
            ]
        );
    }

    #[test]
    fn test_dirichlet_bc() {
        let u1 = CartesianMesh::new(vec![0.0], vec![10.0], vec![5], 1);
        let mut cdf = CartesianDataFrame::new_from(&u1, 1, 2);
        cdf.fill_ic(|x,_,_| x + 1.0);
        cdf.fill_bc(BCType::Dirichlet(0.0), BCType::Dirichlet(1.0));
        assert_eq!(
            cdf.data,
            vec![-6.0, -2.0, 2.0, 4.0, 6.0, 8.0, 10.0, -8.0, -26.0]
        );
    }

    #[test]
    fn test_neumann_bc() {
        let u1 = CartesianMesh::new(vec![0.0], vec![10.0], vec![5], 1);
        let mut cdf = CartesianDataFrame::new_from(&u1, 1, 2);
        cdf.fill_ic(|x, _,_| x + 1.0);
        cdf.fill_bc(BCType::Neumann, BCType::Neumann);
        assert_eq!(
            cdf.data,
            vec![4.0, 2.0, 2.0, 4.0, 6.0, 8.0, 10.0, 10.0, 8.0]
        );
    }

    #[test]
    fn test_prescribed_bc() {
        let u1 = CartesianMesh::new(vec![0.0], vec![10.0], vec![5], 1);
        let mut cdf = CartesianDataFrame::new_from(&u1, 1, 2);
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

