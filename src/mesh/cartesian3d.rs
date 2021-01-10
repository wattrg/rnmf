// Todo:
//     - boundary conditions for any dimension
//     - tests for all dimensions (all tests except boundary conditions have tests for at least
//       1 and 2 dimensions)
//     - labels for each cell, so that cells may be marked as a certain type


use crate::boundary_conditions::*;
use std::rc::Rc;
use crate::Real;
//use super::*;


/// Structure containing data to define a CartesianMesh
#[derive(Debug)]
pub struct CartesianMesh3D {
    /// The position of the nodes in the mesh
    pub node_pos: Vec<Real>,

    /// The lower corner of the mesh
    pub lo: Vec<Real>,

    /// The upper corner of the mesh
    pub hi: Vec<Real>,

    /// The number of cells in nodes in each direction
    pub n : Vec<usize>,

    /// The distance between each node in each
    pub dx: Vec<Real>,

    /// The number of spatial dimensions of the mesh
    pub dim: usize,
    
    /// The number of nodes in the mesh
    n_nodes: usize
}

/// provides functionality to assist with indexing of data. This is separate to the index operator
trait Indexing{
    /// Returns the non-flat index of the item stored at a particular flattened location
    fn get_ijk(&self, p: usize ) -> Option<(isize, isize, isize, usize)>;

    /// Returns a borrow to the value stored at a particular non-flat index
    fn index(&self, i: isize, j: isize, k: isize, n: usize) -> &Real;
}


impl CartesianMesh3D {
    /// Generate new CartesianMesh from the lo and high corner, 
    /// and the number of nodes
    pub fn new(lo: Vec<Real>, hi: Vec<Real>, n:Vec<usize>) -> Rc<CartesianMesh3D>
    {
        let mut dx = vec![0.0; 3];
        for i_dim in 0..3 as usize{
            dx[i_dim] = (hi[i_dim] - lo[i_dim])/(n[i_dim] as Real);
        }

        // calculate the capacity of the vector required to store all the node positions
        let n_nodes = n[0]*n[1]*n[2]*3;

        // allocate memory to the node positions vector
        let node_pos = Vec::with_capacity(n_nodes);

        // allocate memory to the mesh
        let mut cm = CartesianMesh3D
        {
            lo,
            hi,
            n,
            dx,
            node_pos,
            dim: 3,
            n_nodes,
        };

        // calculate the positions of the nodes
        // this will be able to be done better by creating an enumerating iterator
        // which will give (i,j,k,n) to begin with
        cm.node_pos = (0..n_nodes).map(
            |p: usize| -> Real {
                let (i,j,k,n) = cm.get_ijk(p).unwrap();
                match n{
                    0 => {cm.lo[n] + ((i as Real) + 0.5) * cm.dx[n]}
                    1 => {cm.lo[n] + ((j as Real) + 0.5) * cm.dx[n]}
                    2 => {cm.lo[n] + ((k as Real) + 0.5) * cm.dx[n]}
                    _ => {panic!("n cannot be {} in 3D",n)}
                }
                
            }).collect();
        Rc::new(cm)
    }
}

/// function to convert a multi-dimensional coordinate into a single dimension coordinate
fn get_flat_index(i: isize, j: isize, k: isize, n: usize, 
            dimension: &[usize], ng: usize) -> usize {

    let mut p = i + ng as isize;
    p += dimension[0] as isize * (j + ng as isize) 
                + (dimension[0] * dimension[1] * n) as isize;
    
    p += (dimension[0] * dimension[1]) as isize * (k + ng as isize) 
                    + (dimension[2] * 3 * dimension[0] * n) as isize;

    p as usize
}

/// Converts (3+1)D index into a 1D (flattened) index
fn get_ijk_from_p(p: usize, dimension: &[usize], n_nodes: usize, ng: usize) 
                                                        -> Option<(isize, isize, isize, usize)>{
    if p >= n_nodes {
        Option::None
    }

    else{
        let i = (p % dimension[0]) as isize;
        let j: isize;
        let k: isize;
        let d: isize;

        j = (p as isize - i)/dimension[0] as isize % dimension[1] as isize;

        k = (p as isize - j * dimension[1] as isize - i)
                            /((dimension[0] * dimension[1]) as isize) % dimension[2] as isize;
        d = (p as isize - k * (dimension[1] * dimension[0]) as isize - j * dimension[0] as isize - i)
                            /((dimension[0] * dimension[1] * dimension[2]) as isize);
        Some((i - ng as isize, j - ng as isize, k - ng as isize, d as usize))
    }
}

impl Indexing for CartesianMesh3D
{
    /// Retrieves the element at (i,j,k,n)
    fn index(&self, i: isize, j: isize, k: isize, n: usize) -> & Real
    {
        let p = get_flat_index(i, j, k, n, &self.n, 0);
        let np = self.node_pos.get(p as usize);
        match np{
            Some(np) => {np},
            None => panic!("Index ({}, {}, {}, {}) out of range of Cartesian mesh with size {:?}", 
                                i,j,k,n, self.n),
        }

    }

    /// Returns the un-flattened index
    fn get_ijk(& self, p: usize) -> Option<(isize, isize, isize, usize)>{
        get_ijk_from_p(p, &self.n, self.n_nodes, 0)
    }
}


/// Immutable Iterator struct for CartesianMesh
pub struct CartesianMeshIter <'a> {
    mesh: &'a CartesianMesh3D,
    current_indx: usize,
}

impl <'a> Iterator for CartesianMeshIter<'a>{
    // return a vector of references instead of slice because the length of the returned object
    // depends on the number of dimensions, and thus is not known at compile time, so slices
    // won't work
    type Item = Vec<&'a Real>;
    
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

impl <'a> IntoIterator for &'a CartesianMesh3D {
    type Item = Vec<&'a Real>;
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
#[derive(Debug, Clone)]
pub struct CartesianDataFrame3D{
    /// The data is stored here
    pub data: Vec<Real>,

    /// The number of ghost nodes, added to each side of the underlying `CartesianMesh`
    pub n_ghost: usize,

    n_grown: Vec<usize>,

    /// Reference to the underlying `CartesianMesh`
    pub underlying_mesh: Rc<CartesianMesh3D>,

    /// The number of components to be stored in the data frame. For example, storing
    /// a 3D velocity would result in `n_comp = 3`
    pub n_comp: usize,

    /// The total number of individual pieces of information needed to be stored
    n_nodes: usize,
}


/// data structure to store data on CartesianMesh
impl CartesianDataFrame3D {
    /// generate new `CartesianDataFrame` from a given mesh, adding a given
    /// number of ghost nodes
    pub fn new_from(m: & Rc<CartesianMesh3D>, 
                    n_comp: usize, n_ghost: usize) -> CartesianDataFrame3D
    {
        let mut n_nodes = m.n[0] + 2*n_ghost;
        if m.dim >= 2 {
            n_nodes *= m.n[1] + 2*n_ghost;

            if m.dim == 3 {
                n_nodes *= m.n[2] + 2*n_ghost;
            }
        }
        n_nodes *= n_comp;


        CartesianDataFrame3D
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
    pub fn fill_ic (&mut self, ic: impl Fn(Real, Real, Real, usize)->Real)
    {
        for (pos, val) in self.into_iter().enumerate_pos(){
            let (x,y,z,n) = pos;
            *val = ic(x,y,z,n);
        }

    }

    /// Convert the flat index of a data point into the (3+1)D coordinates
    fn get_ijk(&self, p: usize) -> Option<(isize, isize, isize, usize)> {
        get_ijk_from_p(p, &self.n_grown, self.n_nodes, self.n_ghost)
    }
}

impl core::ops::IndexMut<(isize, isize, isize, usize)> for CartesianDataFrame3D {
    /// Exclusively borrows element at (i,j,k,n). The valid cells are indexed from zero 
    /// and ghost cells at the lower side of the domain are indexed with negative numbers
    fn index_mut(&mut self, indx: (isize, isize,isize,usize)) -> &mut Real {
        let (i,j,k,n) = indx;
        let p = get_flat_index(i,j,k,n,&self.n_grown,self.n_ghost);
        &mut self.data[p]
    }
}


impl core::ops::Index<(isize, isize, isize, usize)> for CartesianDataFrame3D{
    type Output = Real;

    #[allow(dead_code)] 
    /// Borrows the element at (i,j,k,n). The valid cells are indexed from zero
    /// and ghost cells at the lower side of the domain are indexed with negative
    /// numbers.
    fn index(&self, indx: (isize, isize, isize, usize) ) -> &Real {
        let (i,j,k,n) = indx;
        let p = get_flat_index(i, j, k, n, &self.n_grown, self.n_ghost);

        &self.data[p]
    }
}


impl <'a> BoundaryCondition for CartesianDataFrame3D{
    /// Fill the ghost nodes of the CartesianDataFrame based on BCType
    fn fill_bc (&mut self, bc: &BCs) {
        for i_comp in 0..self.n_comp{
            for i_dim in 0..self.underlying_mesh.dim{
                let bc_lo = &bc.bcs[i_comp].lo[i_dim];
                let bc_hi = &bc.bcs[i_comp].hi[i_dim];
                // low boundary condition
                match bc_lo {
                    BCType::Prescribed(values) => {
                        panic!("Prescribed BC in 3D not supported yet");
                    }
                    BCType::Neumann(gradient) => {
                        panic!("Neumann BC in 3D not supported yet");
                    }
                    BCType::Dirichlet(value) => {
                        panic!("Dirichlet BC for 3D not yet supported");
                    }
                }

                // high boundary condition
                match bc_hi {
                    BCType::Prescribed(values) => {
                        panic!("Prescribed BC in 3D not supported yet");
                    }
                    BCType::Neumann(gradient) => {
                        panic!("Neumann BC in 3D not supported yet");
                    }
                    BCType::Dirichlet(value) => {
                        panic!("Dirichlet BC for 3D not yet supported");
                    }
                }
            }
        }
    }

    /// Determines if the cell at (i,j,k) is a valid cell (returns true) or a ghost
    /// cell (returns false)
    fn ijk_is_valid_cell(&self, i: isize, j: isize, k: isize) -> bool {
        let above_zero = i >= 0 && j >= 0 && k >= 0;

        let mut below_max = i < self.underlying_mesh.n[0] as isize;
        below_max = below_max && j < self.underlying_mesh.n[1]  as isize;
        below_max = below_max && k < self.underlying_mesh.n[2] as isize;
        above_zero && below_max
    }

    fn p_is_valid_cell(&self, p: usize) -> bool {
        let (i,j,k,_) = self.get_ijk(p).unwrap();
        self.ijk_is_valid_cell(i, j, k)
    }
}



/// mutable iterator for `CartesianDataFrame`.
pub struct CartesianDataFrameIter<'a> {
    df: &'a mut CartesianDataFrame3D,
    current_indx: usize,
}



impl <'a> Iterator for CartesianDataFrameIter<'a> {
    type Item = &'a mut Real;

    // this is a safe function wrapping unsafe code. Rust cannot guarantee that it is safe, but
    // in practice it can be, and I'm pretty sure that it should be safe for everything we will 
    // ever want to do
    fn next(&mut self) -> Option<Self::Item> {
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

impl <'a> IntoIterator for &'a mut CartesianDataFrame3D{
    type Item = &'a mut Real;
    type IntoIter = CartesianDataFrameIter<'a>;

    fn into_iter(self) -> Self::IntoIter{
        Self::IntoIter {
            current_indx: 0,
            df: self,
        }
    }
}

/// Struct to iterate over a data frame, with the current index and the data itself
pub struct EnumerateIndex<'a>{
    iter: CartesianDataFrameIter<'a>,
}

/// Turns `CartesianDataFrameIter` into an iterator which enumerates over the current index
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
    type Item = ((isize, isize, isize, usize), &'a mut Real);

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

/// Structure to iterate over data, giving both the position of the data, and the data itself
pub struct EnumeratePos<'a>{
    iter: CartesianDataFrameIter<'a>,
}

/// Turns `CartesianDataFrameIter` into an iterator which enumerates the position of the data
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
    type Item = ((Real, Real, Real, usize), &'a mut Real);

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

impl std::ops::Mul<CartesianDataFrame3D> for Real {
    type Output = CartesianDataFrame3D;

    fn mul(self, rhs: CartesianDataFrame3D) -> Self::Output {
        let mut result = CartesianDataFrame3D::new_from(&rhs.underlying_mesh, rhs.n_comp, rhs.n_ghost);
        for (i, vals) in result.data.iter_mut().enumerate(){
            *vals = self * rhs.data[i];
        }
        result
    }
}
impl std::ops::Mul<&CartesianDataFrame3D> for Real {
    type Output = CartesianDataFrame3D;

    fn mul(self, rhs: &CartesianDataFrame3D) -> Self::Output {
        let mut result = CartesianDataFrame3D::new_from(&rhs.underlying_mesh, rhs.n_comp, rhs.n_ghost);
        for (i, vals) in result.data.iter_mut().enumerate(){
            *vals = self * rhs.data[i];
        }
        result
    }
}

impl std::ops::Add<&CartesianDataFrame3D> for &CartesianDataFrame3D {
    type Output = CartesianDataFrame3D;

    fn add(self, rhs: &CartesianDataFrame3D) -> Self::Output {
        let mut sum = CartesianDataFrame3D::new_from(&rhs.underlying_mesh, rhs.n_comp, rhs.n_ghost);
        for (i, vals) in sum.data.iter_mut().enumerate() {
            *vals = self.data[i] + rhs.data[i];
        }
        sum
    }
}

impl std::ops::Mul<&CartesianDataFrame3D> for &CartesianDataFrame3D {
    type Output = CartesianDataFrame3D;

    fn mul(self, rhs: &CartesianDataFrame3D) -> Self::Output {
        let mut result = CartesianDataFrame3D::new_from(&rhs.underlying_mesh, rhs.n_comp, rhs.n_ghost);
        for (i, vals) in result.data.iter_mut().enumerate(){
            *vals = self.data[i] * rhs.data[i];
        }
        result
    }
}




//impl DataFrame for CartesianDataFrame {}


/// mutable iterator for `CartesianDataFrame`.
pub struct CartesianDataFrameGhostIter<'a> {
    df: &'a mut CartesianDataFrame3D,
    current_indx: usize,
}



impl <'a> Iterator for CartesianDataFrameGhostIter<'a> {
    type Item = &'a mut Real;

    // this is a safe function wrapping unsafe code. Rust cannot guarantee that it is safe, but
    // in practice it can be, and I'm pretty sure that it should be safe for everything we will 
    // ever want to do
    fn next(&mut self) -> Option<Self::Item> {
        // progress the current index to skip ghost cells
        while self.current_indx <= self.df.n_nodes - self.df.n_comp &&
                                   self.df.p_is_valid_cell(self.current_indx) {
            self.current_indx += 1;
        }

        // check if the next item exists
        if self.current_indx <= self.df.n_nodes-self.df.n_comp 
                   && !self.df.p_is_valid_cell(self.current_indx){
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

pub trait GhostIterator <'a> {
    fn ghost(self) -> CartesianDataFrameGhostIter<'a>;
}

impl <'a> GhostIterator <'a> for CartesianDataFrameIter<'a> {
    fn ghost(self) -> CartesianDataFrameGhostIter<'a> {
        CartesianDataFrameGhostIter{
            current_indx: 0,
            df: self.df,
        }   
    }
}


/// Struct to iterate over a data frame, with the current index and the data itself
pub struct EnumerateGhostIndex<'a>{
    iter: CartesianDataFrameGhostIter<'a>,
}

/// Turns `CartesianDataFrameIter` into an iterator which enumerates over the current index
pub trait GhostIndexEnumerable <'a> {
    fn enumerate_index(self) -> EnumerateGhostIndex<'a>;
}

impl <'a> GhostIndexEnumerable <'a> for CartesianDataFrameGhostIter<'a> {
    fn enumerate_index(self) -> EnumerateGhostIndex<'a> {
        EnumerateGhostIndex{
            iter: self,
        }
    }
}

impl <'a> Iterator for EnumerateGhostIndex<'a> {
    type Item = ((isize, isize, isize, usize), &'a mut Real);

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



#[cfg(test)]
mod tests{
    use super::*;






    #[test]
    fn indexing () {

        let m3 = CartesianMesh3D::new(vec![0.0,0.0,0.0], vec![6.0,6.0,6.0], vec![3,3,3]);
        assert_eq!(m3.get_ijk(3).unwrap(), (0,1,0,0));
        assert_eq!(m3.get_ijk(22).unwrap(), (1,1,2,0));
        assert_eq!(m3.get_ijk(81), None);
        assert_eq!(m3.index(2,1,0,1), &3.0);

    }

    #[test]
    fn node_pos () {
        let m3 = CartesianMesh3D::new(vec![0.0,0.0,0.0], vec![6.0,6.0,6.0], vec![3,3,3]);
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
}