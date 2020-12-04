// Todo:
//     - boundary conditions for any dimension
//     - tests for all dimensions (all tests except boundary conditions have tests for at least
//       1 and 2 dimensions)


use crate::boundary_conditions::{BCType, BoundaryCondition};
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
    fn get_ijk(&self, p: usize ) -> Option<(isize, isize, isize, usize)>;

    /// Returns a borrow to the value stored at a particular non-flat index
    fn index(&self, i: isize, j: isize, k: isize, n: usize) -> &f64;
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

/// function to convert a multi-dimensional coordinate into a single dimension coordinate
fn get_flat_index(i: isize, j: isize, k: isize, n: usize, 
            dimension: &[usize], dim: usize, ng: usize) -> usize {

    let mut p = i + ng as isize;
    if dim >= 2 {
        p += dimension[0] as isize * (j + ng as isize) 
                + (dimension[0] * dimension[1] * n) as isize;
    
        if dim == 3 {
            p += (dimension[0] * dimension[1]) as isize * (k + ng as isize) 
                    + (dimension[2] * dim * dimension[0] * n) as isize;
        }
    }

    p as usize
}

fn get_ijk_from_p(p: usize, dimension: &[usize], n_nodes: usize, dim: usize, ng: usize) 
                                                        -> Option<(isize, isize, isize, usize)>{
    if p >= n_nodes {
        Option::None
    }

    else{
        let i = (p % dimension[0]) as isize;
        let mut j = 0;
        let k: isize;
        let d: isize;

        if dim >= 2 { // either 2D or 3D
            j = (p as isize - i)/dimension[0] as isize % dimension[1] as isize;
        }

        match dim{
            1 => {
                d = (p as isize - i)/dimension[0] as isize;
                Some((i - ng as isize, 0, 0, d as usize))
            }
            2 => {
                d = (p as isize - j * dimension[1] as isize - i)
                            /dimension[0] as isize % dim as isize;
                Some((i - ng as isize, j - ng as isize, 0, d as usize))
            }
            3 => {
                k = (p as isize - j * dimension[1] as isize - i)
                            /((dimension[0] * dimension[1]) as isize) % dimension[2] as isize;
                d = (p as isize - k * (dimension[1] * dimension[0]) as isize - j * dimension[0] as isize - i)
                            /((dimension[0] * dimension[1] * dimension[2]) as isize);
                Some((i - ng as isize, j - ng as isize, k - ng as isize, d as usize))
            }
            _ => {
                panic!("{}D not supported!", dim);
            }
        }
    }
}

impl Indexing for CartesianMesh
{
    /// Retrieves the element at (i,j,k,n)
    #[allow(dead_code)]
    fn index(&self, i: isize, j: isize, k: isize, n: usize) -> & f64
    {
        let p = get_flat_index(i, j, k, n, &self.n, self.dim, 0);
        let np = self.node_pos.get(p as usize);
        match np{
            Some(np) => {np},
            None => panic!("Index ({}, {}, {}, {}) out of range of Cartesian mesh with size {:?}", 
                                i,j,k,n, self.n),
        }

    }

    /// Returns the un-flattened index
    fn get_ijk(& self, p: usize) -> Option<(isize, isize, isize, usize)>{
        get_ijk_from_p(p, &self.n, self.n_nodes, self.dim, 0)
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
        if self.current_indx >= self.mesh.node_pos.len(){
            None
        }
        // If the next position exists, return a vector of size self.dim containing references
        // to the position
        else {
            let mut next_pos = Vec::with_capacity(self.mesh.dim);
            let (i,j,k,n) = self.mesh.get_ijk(self.current_indx).unwrap();
            for i_dim in 0..(self.mesh.dim) {
                next_pos.push(self.mesh.index(i, j, k, n + i_dim))
            }
            self.current_indx += 1;
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
        let mut n_nodes = m.n[0] + 2*n_ghost;
        if m.dim >= 2 {
            n_nodes *= m.n[1] + 2*n_ghost;

            if m.dim == 3 {
                n_nodes *= m.n[2] + 2*n_ghost;
            }
        }
        n_nodes *= n_comp;


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
    pub fn fill_ic (&mut self, _ic: fn(f64, f64, f64, usize)->f64)
    {
        for (pos, val) in self.into_iter().enumerate_pos(){
            let (x,y,z,n) = pos;
            *val = _ic(x,y,z,n);
        }

    }

    /// Convert the flat index of a data point into the (3+1)D coordinates
    fn get_ijk(&self, p: usize) -> Option<(isize, isize, isize, usize)> {
        get_ijk_from_p(p, &self.n_grown, self.n_nodes, self.underlying_mesh.dim, self.n_ghost)
    }
}

impl core::ops::IndexMut<(isize, isize, isize, usize)> for CartesianDataFrame {
    /// Exclusively borrows element at (i,j,k,n). The valid cells are indexed from zero 
    /// and ghost cells at the lower side of the domain are indexed with negative numbers
    fn index_mut(&mut self, indx: (isize, isize,isize,usize)) -> &mut f64 {
        let (i,j,k,n) = indx;
        let p = get_flat_index(i,j,k,n,&self.n_grown,self.underlying_mesh.dim,self.n_ghost);
        &mut self.data[p]
    }
}


impl core::ops::Index<(isize, isize, isize, usize)> for CartesianDataFrame{
    type Output = f64;

    #[allow(dead_code)] 
    /// Borrows the element at (i,j,k,n). The valid cells are indexed from zero
    /// and ghost cells at the lower side of the domain are indexed with negative
    /// numbers.
    fn index(&self, indx: (isize, isize, isize, usize) ) -> &f64 {
        let (i,j,k,n) = indx;
        let p = get_flat_index(i, j, k, n, &self.n_grown, self.underlying_mesh.dim, self.n_ghost);

        &self.data[p]
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

    fn ijk_is_valid_cell(&self, i: isize, j: isize, k: isize) -> bool {
        let above_zero = i >= 0 && j >= 0 && k >= 0;

        let mut below_max = i < self.underlying_mesh.n[0] as isize;
        if self.underlying_mesh.dim >= 2 {
            below_max = below_max && j < self.underlying_mesh.n[1]  as isize;
            if self.underlying_mesh.dim == 3 {
                below_max = below_max && k < self.underlying_mesh.n[2] as isize;
            }
        }
        above_zero && below_max
    }

    fn p_is_valid_cell(&self, p: usize) -> bool {
        let (i,j,k,_) = self.get_ijk(p).unwrap();
        self.ijk_is_valid_cell(i, j, k)
    }
}



/// mutable iterator for `CartesianDataFrame`.
pub struct CartesianDataFrameIter<'a> {
    df: &'a mut CartesianDataFrame,
    current_indx: usize,
}



impl <'a> Iterator for CartesianDataFrameIter<'a> {
    type Item = &'a mut f64;

    // this is a safe function wrapping unsafe code. Rust cannot guarantee that it is safe, but
    // in practice it can be, and I'm pretty sure that it should be safe for everything we will 
    // ever want to do
    fn next(&mut self) -> Option<Self::Item> {
        //println!("Calling next on base iterator");
        // progress the current index to skip ghost cells
        while self.current_indx <= self.df.n_nodes - self.df.n_comp &&
                                  !self.df.p_is_valid_cell(self.current_indx) {
            self.current_indx += 1;
        }

        // check if the next item exists
        if self.current_indx <= self.df.n_nodes-self.df.n_comp 
                   && self.df.p_is_valid_cell(self.current_indx){
            // access and return next item 
            let ptr = self.df.data.as_mut_ptr();
            let next_data: Option<Self::Item>;
            unsafe{
                next_data = ptr.add(self.current_indx).as_mut();
            }
            self.current_indx += 1;
            next_data
        }
        // if the next item doesn't exist, return None
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

pub struct EnumerateIndex<'a>{
    iter: CartesianDataFrameIter<'a>,
}

pub trait IndexEnumerable <'a> {
    fn enumerate_index(self) -> EnumerateIndex<'a>;
}

impl <'a> IndexEnumerable <'a> for CartesianDataFrameIter<'a> {
    fn enumerate_index(self) -> EnumerateIndex<'a> {
        EnumerateIndex{
            iter: self,
        }
    }
}

impl <'a> Iterator for EnumerateIndex<'a> {
    type Item = ((isize, isize, isize, usize), &'a mut f64);

    fn next(&mut self) -> Option<Self::Item>{
        let next_item = self.iter.next();
        match next_item{
            Some(nv) => {
                Some((self.iter.df.get_ijk(self.iter.current_indx-1).unwrap(), nv))
            }
            None => {
                None
            }
        }
        
    }
}

pub struct EnumeratePos<'a>{
    iter: CartesianDataFrameIter<'a>,
}
pub trait PosEnumerable <'a> {
    fn enumerate_pos(self) -> EnumeratePos<'a>;
}

impl <'a> PosEnumerable <'a> for CartesianDataFrameIter<'a> {
    fn enumerate_pos(self) -> EnumeratePos<'a> {
        EnumeratePos{
            iter: self,
        }
    }
}

impl <'a> Iterator for EnumeratePos<'a> {
    type Item = ((f64, f64, f64, usize), &'a mut f64);

    fn next(&mut self) -> Option<Self::Item>{
        let next_item = self.iter.next();
        match next_item{
            Some(nv) => {
                let (i,j,k,n) = self.iter.df.get_ijk(self.iter.current_indx-1).unwrap();
                let pos_x = *self.iter.df.underlying_mesh.index(i, j, k, 0);
                let mut pos_y = 0.0;
                let mut pos_z = 0.0;
                if self.iter.df.underlying_mesh.dim >= 2 {
                    pos_y = *self.iter.df.underlying_mesh.index(i, j, k, 1);
                    if self.iter.df.underlying_mesh.dim == 3 {
                        pos_z = *self.iter.df.underlying_mesh.index(i, j, k, 2);
                    }
                }
                Some(((pos_x, pos_y, pos_z, n), nv))
            }
            None => {
                None
            }
        }
        
    }
}



#[cfg(test)]
mod tests{
    use super::*;

    #[test]
    fn is_valid_cell () {
        let m1 = CartesianMesh::new(vec![0.0], vec![6.0], vec![3], 1);
        let df = CartesianDataFrame::new_from(&m1, 1, 1);

        assert_eq!(df.p_is_valid_cell(0), false);
        assert_eq!(df.p_is_valid_cell(1), true);
        assert_eq!(df.p_is_valid_cell(2), true);
        assert_eq!(df.p_is_valid_cell(3), true);
        assert_eq!(df.p_is_valid_cell(4), false);
    }

    #[test]
    fn enumerate_pos () {
        let m1 = CartesianMesh::new(vec![0.0], vec![6.0], vec![3], 1);
        let mut df = CartesianDataFrame::new_from(&m1, 1, 1);
        

        let mut df_iter_pos = df.into_iter().enumerate_pos();

        assert_eq!(df_iter_pos.next(), Some(((1.0, 0.0, 0.0, 0), &mut 0.0)));
        assert_eq!(df_iter_pos.next(), Some(((3.0, 0.0, 0.0, 0), &mut 0.0)));
        

    }

    #[test]
    // currently only tests the iterator on 1D data
    fn data_frame_iterator() {
        // 1D
        {
            let m1 = CartesianMesh::new(vec![0.0], vec![6.0], vec![3], 1);
            let mut df = CartesianDataFrame::new_from(&m1, 1, 1);
            df.fill_ic(|x,_,_,_|  x + 1.0);

            // test that the iterators gives the correct values
            let mut df_iter = df.into_iter();
            assert_eq!(df_iter.next(), Some(&mut 2.0));
            assert_eq!(df_iter.next(), Some(&mut 4.0));
            assert_eq!(df_iter.next(), Some(&mut 6.0));
            assert_eq!(df_iter.next(), Option::None);
            
            // test to make sure mutating the data works
            for data in &mut df{
                *data += 1.0;
            }
            assert_eq!(df.data, vec![0.0, 3.0, 5.0, 7.0, 0.0]);
        }
        
        // 2D
        {
            let m2 = CartesianMesh::new(vec![0.0, 0.0], vec![6.0, 6.0], vec![3, 3], 2);
            let mut df2 = CartesianDataFrame::new_from(&m2, 1, 1);
            df2.fill_ic(|x,y,_,_| x + y);
            


            let mut df2_iter = df2.into_iter();
            assert_eq!(df2_iter.next(), Some(&mut 2.0));
            assert_eq!(df2_iter.next(), Some(&mut 4.0));
            assert_eq!(df2_iter.next(), Some(&mut 6.0));
            assert_eq!(df2_iter.next(), Some(&mut 4.0));
            assert_eq!(df2_iter.next(), Some(&mut 6.0));
            assert_eq!(df2_iter.next(), Some(&mut 8.0));
            assert_eq!(df2_iter.next(), Some(&mut 6.0));
            assert_eq!(df2_iter.next(), Some(&mut 8.0));
            assert_eq!(df2_iter.next(), Some(&mut 10.0));
            assert_eq!(df2_iter.next(), None);
        }
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
        assert_eq!(m1.index(3, 0, 0, 0), &7.0);
        let mut cdf1 = CartesianDataFrame::new_from(&m1, 1, 2);
        cdf1.fill_ic(|x,_,_,_|x + 1.0);

        let m2 = CartesianMesh::new(vec![0.0,0.0], vec![6.0,6.0], vec![3,3], 2);
        assert_eq!(m2.get_ijk(16).unwrap(), (1,2,0,1));
        assert_eq!(m2.get_ijk(40), Option::None);
        assert_eq!(m2.index(2,1,0,1), &3.0);

        let m3 = CartesianMesh::new(vec![0.0,0.0,0.0], vec![6.0,6.0,6.0], vec![3,3,3], 3);
        assert_eq!(m3.get_ijk(3).unwrap(), (0,1,0,0));
        assert_eq!(m3.get_ijk(22).unwrap(), (1,1,2,0));
        assert_eq!(m3.get_ijk(81), None);
        assert_eq!(m3.index(2,1,0,1), &3.0);
    }

    #[test]
    fn node_pos () {
        let m1 = CartesianMesh::new(vec![0.0], vec![10.0], vec![5], 1);
        assert_eq!(m1.node_pos, vec![1.0, 3.0, 5.0, 7.0, 9.0]);

        let m2 = CartesianMesh::new(vec![0.0, 0.0], vec![10.0, 10.0], vec![5,5], 2);
        assert_eq!(
            m2.node_pos,
            vec![
                // x values
                1.0, 3.0, 5.0, 7.0, 9.0, 
                1.0, 3.0, 5.0, 7.0, 9.0,
                1.0, 3.0, 5.0, 7.0, 9.0,
                1.0, 3.0, 5.0, 7.0, 9.0,
                1.0, 3.0, 5.0, 7.0, 9.0,
                // y values
                1.0, 1.0, 1.0, 1.0, 1.0,
                3.0, 3.0, 3.0, 3.0, 3.0, 
                5.0, 5.0, 5.0, 5.0, 5.0,
                7.0, 7.0, 7.0, 7.0, 7.0, 
                9.0, 9.0, 9.0, 9.0, 9.0
            ]
        );

        let m3 = CartesianMesh::new(vec![0.0,0.0,0.0], vec![6.0,6.0,6.0], vec![3,3,3], 3);
        assert_eq!(
            m3.node_pos,
            vec![
                // x values
                1.0, 3.0, 5.0, 
                1.0, 3.0, 5.0, 
                1.0, 3.0, 5.0,

                1.0, 3.0, 5.0, 
                1.0, 3.0, 5.0, 
                1.0, 3.0, 5.0,

                1.0, 3.0, 5.0, 
                1.0, 3.0, 5.0, 
                1.0, 3.0, 5.0,
                
                // y values
                1.0, 1.0, 1.0, 
                3.0, 3.0, 3.0, 
                5.0, 5.0, 5.0,

                1.0, 1.0, 1.0, 
                3.0, 3.0, 3.0, 
                5.0, 5.0, 5.0,

                1.0, 1.0, 1.0, 
                3.0, 3.0, 3.0, 
                5.0, 5.0, 5.0,
                
                // z values
                1.0, 1.0, 1.0, 
                1.0, 1.0, 1.0, 
                1.0, 1.0, 1.0,

                3.0, 3.0, 3.0, 
                3.0, 3.0, 3.0, 
                3.0, 3.0, 3.0,

                5.0, 5.0, 5.0, 
                5.0, 5.0, 5.0, 
                5.0, 5.0, 5.0,
            ]
        );
    }

    #[test]
    fn test_dirichlet_bc() {
        let u1 = CartesianMesh::new(vec![0.0], vec![10.0], vec![5], 1);
        let mut cdf = CartesianDataFrame::new_from(&u1, 1, 2);
        cdf.fill_ic(|x,_,_,_| x + 1.0);
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
        cdf.fill_ic(|x, _,_,_| x + 1.0);
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
        cdf.fill_ic(|x, _,_,_| x + 1.0);
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

