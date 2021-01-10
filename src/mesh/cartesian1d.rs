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
pub struct CartesianMesh1D {
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

    /// The number of nodes in the mesh
    n_nodes: usize,

    pub dim: usize,
}

/// provides functionality to assist with indexing of data. This is separate to the index operator
trait Indexing{
    /// Returns the non-flat index of the item stored at a particular flattened location
    fn get_ijk(&self, p: usize ) -> Option<(isize, usize)>;

    /// Returns a borrow to the value stored at a particular non-flat index
    fn index(&self, i: isize, n: usize) -> &Real;
}


impl CartesianMesh1D {
    /// Generate new CartesianMesh from the lo and high corner, 
    /// and the number of nodes
    pub fn new(lo: Vec<Real>, hi: Vec<Real>, n:Vec<usize>) -> Rc<CartesianMesh1D>
    {
        let mut dx = vec![0.0; 1];
        dx[0] = (hi[0] - lo[0])/(n[0] as Real);

        // calculate the capacity of the vector required to store all the node positions
        let n_nodes = n[0];

        // allocate memory to the node positions vector
        let node_pos = Vec::with_capacity(n_nodes);

        // allocate memory to the mesh
        let mut cm = CartesianMesh1D
        {
            lo,
            hi,
            n,
            dx,
            node_pos,
            n_nodes,
            dim: 1,
        };

        // calculate the positions of the nodes
        // this will be able to be done better by creating an enumerating iterator
        // which will give (i,j,k,n) to begin with
        cm.node_pos = (0..n_nodes).map(
            |p: usize| -> Real {
                let (i,n) = cm.get_ijk(p).unwrap();
                cm.lo[n] + ((i as Real) + 0.5) * cm.dx[n]
                
            }).collect();
        Rc::new(cm)
    }
}

/// function to convert a multi-dimensional coordinate into a single dimension coordinate
fn get_flat_index(i: isize, n: usize, dimension: &[usize], ng: usize) -> usize {

    i as usize + ng + dimension[0]*n

}

/// Converts (3+1)D index into a 1D (flattened) index
fn get_ijk_from_p(p: usize, dimension: &[usize], n_nodes: usize, ng: usize) 
                                                        -> Option<(isize, usize)>{
    if p >= n_nodes {
        Option::None
    }

    else{
        let i = (p % dimension[0]) as isize;
        let d: isize;

        d = (p as isize - i)/dimension[0] as isize;
        Some((i - ng as isize, d as usize))
    }
}

impl Indexing for CartesianMesh1D
{
    /// Retrieves the element at (i,j,k,n)
    fn index(&self, i: isize, n: usize) -> & Real
    {
        let p = get_flat_index(i, n, &self.n, 0);
        let np = self.node_pos.get(p as usize);
        match np{
            Some(np) => {np},
            None => panic!("Index ({}, {}) out of range of Cartesian mesh with size {:?}", 
                                i,n, self.n),
        }

    }

    /// Returns the un-flattened index
    fn get_ijk(& self, p: usize) -> Option<(isize, usize)>{
        get_ijk_from_p(p, &self.n, self.n_nodes, 0)
    }
}


/// Immutable Iterator struct for CartesianMesh
pub struct CartesianMeshIter <'a> {
    mesh: &'a CartesianMesh1D,
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
            let (i,n) = self.mesh.get_ijk(self.current_indx).unwrap();
            for i_dim in 0..(self.mesh.dim) {
                next_pos.push(self.mesh.index(i, n + i_dim))
            }
            self.current_indx += 1;
            Some(next_pos)
        }
    }
}

impl <'a> IntoIterator for &'a CartesianMesh1D {
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
pub struct CartesianDataFrame1D{
    /// The data is stored here
    pub data: Vec<Real>,

    /// The number of ghost nodes, added to each side of the underlying `CartesianMesh`
    pub n_ghost: usize,

    n_grown: Vec<usize>,

    /// Reference to the underlying `CartesianMesh`
    pub underlying_mesh: Rc<CartesianMesh1D>,

    /// The number of components to be stored in the data frame. For example, storing
    /// a 3D velocity would result in `n_comp = 3`
    pub n_comp: usize,

    /// The total number of individual pieces of information needed to be stored
    n_nodes: usize,
}


/// data structure to store data on CartesianMesh
impl CartesianDataFrame1D {
    /// generate new `CartesianDataFrame` from a given mesh, adding a given
    /// number of ghost nodes
    pub fn new_from(m: & Rc<CartesianMesh1D>, 
                    n_comp: usize, n_ghost: usize) -> CartesianDataFrame1D
    {
        let mut n_nodes = m.n[0] + 2*n_ghost;
        n_nodes *= n_comp;


        CartesianDataFrame1D
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
    pub fn fill_ic (&mut self, ic: impl Fn(Real, usize)->Real)
    {
        for (pos, val) in self.into_iter().enumerate_pos(){
            let (x,n) = pos;
            *val = ic(x,n);
        }

    }

    /// Convert the flat index of a data point into the (3+1)D coordinates
    fn get_ijk(&self, p: usize) -> Option<(isize, usize)> {
        get_ijk_from_p(p, &self.n_grown, self.n_nodes, self.n_ghost)
    }
}

impl core::ops::IndexMut<(isize, usize)> for CartesianDataFrame1D {
    /// Exclusively borrows element at (i,j,k,n). The valid cells are indexed from zero 
    /// and ghost cells at the lower side of the domain are indexed with negative numbers
    fn index_mut(&mut self, indx: (isize,usize)) -> &mut Real {
        let (i,n) = indx;
        let p = get_flat_index(i,n,&self.n_grown,self.n_ghost);
        &mut self.data[p]
    }
}


impl core::ops::Index<(isize, usize)> for CartesianDataFrame1D{
    type Output = Real;

    #[allow(dead_code)] 
    /// Borrows the element at (i,j,k,n). The valid cells are indexed from zero
    /// and ghost cells at the lower side of the domain are indexed with negative
    /// numbers.
    fn index(&self, indx: (isize, usize) ) -> &Real {
        let (i,n) = indx;
        let p = get_flat_index(i, n, &self.n_grown, self.n_ghost);

        &self.data[p]
    }
}


pub trait BoundaryCondition1D {
    /// function which fills the boundary conditions
    fn fill_bc(&mut self, bcs: &BCs);

    /// checks if the cell at (i,j,k) contains a valid or a ghost cell. Returns true if valid,
    /// and returns false if ghost
    fn ijk_is_valid_cell(&self, i: isize) -> bool;

    /// check if the cell at p contains a valid or ghost cell. Returns same as ijk_is_valid_cell
    fn p_is_valid_cell(&self, p: usize) -> bool;
}
impl <'a> BoundaryCondition1D for CartesianDataFrame1D{
    /// Fill the ghost nodes of the CartesianDataFrame based on BCType
    fn fill_bc (&mut self, bc: &BCs) {
        for i_comp in 0..self.n_comp{
            let bc_lo = &bc.bcs[i_comp].lo[0];
            let bc_hi = &bc.bcs[i_comp].hi[0];

            // low boundary condition
            match bc_lo {
                BCType::Prescribed(values) => {
                    for (i, &val) in values.iter().rev().enumerate() {
                        self.data[i] = val;
                    }
                        
                }
                BCType::Neumann(gradient) => {
                        for i in 0usize..self.n_ghost as usize {
                            self.data[self.n_ghost - i - 1] = self.data[self.n_ghost + i];
                        }
                    }
                BCType::Dirichlet(value) => {
                    let m: Real = self.data[self.n_ghost as usize] - value;
                    for i in 0usize..self.n_ghost as usize {
                        self.data[self.n_ghost - i - 1] =
                            self.data[self.n_ghost] - 2.0 * (i as Real + 1.0) * m
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
            
                BCType::Neumann(gradient) => {
                    for i in 1usize..self.n_ghost as usize + 1 {
                        self.data[n - (self.n_ghost as usize) + i] =
                            self.data[n - self.n_ghost as usize - i + 1];
                    }
                }
                BCType::Dirichlet(value) => {
                    let m = value - self.data[n - self.n_ghost as usize];
                    for i in 1usize..self.n_ghost as usize + 1 {
                        self.data[n - (self.n_ghost as usize) + i] =
                            self.data[n - self.n_ghost as usize] + 2.0 * m * i as Real;
                    }
                }
            }
        }
    }

    /// Determines if the cell at (i,j,k) is a valid cell (returns true) or a ghost
    /// cell (returns false)
    fn ijk_is_valid_cell(&self, i: isize) -> bool {
        let above_zero = i >= 0;

        let below_max = i < self.underlying_mesh.n[0] as isize;
        above_zero && below_max
    }

    fn p_is_valid_cell(&self, p: usize) -> bool {
        let (i,_) = self.get_ijk(p).unwrap();
        self.ijk_is_valid_cell(i)
    }
}




/// mutable iterator for `CartesianDataFrame`.
pub struct CartesianDataFrameIter<'a> {
    df: &'a mut CartesianDataFrame1D,
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

impl <'a> IntoIterator for &'a mut CartesianDataFrame1D{
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
    type Item = ((isize, usize), &'a mut Real);

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
    type Item = ((Real, usize), &'a mut Real);

    fn next(&mut self) -> Option<Self::Item>{
        let next_item = self.iter.next();
        match next_item{
            Some(nv) => {
                let (i,n) = self.iter.df.get_ijk(self.iter.current_indx-1).unwrap();
                let pos_x = *self.iter.df.underlying_mesh.index(i, 0);
                Some(((pos_x, n), nv))
            }
            None => {
                None
            }
        }
        
    }
}

impl std::ops::Mul<CartesianDataFrame1D> for Real {
    type Output = CartesianDataFrame1D;

    fn mul(self, rhs: CartesianDataFrame1D) -> Self::Output {
        let mut result = CartesianDataFrame1D::new_from(&rhs.underlying_mesh, rhs.n_comp, rhs.n_ghost);
        for (i, vals) in result.data.iter_mut().enumerate(){
            *vals = self * rhs.data[i];
        }
        result
    }
}
impl std::ops::Mul<&CartesianDataFrame1D> for Real {
    type Output = CartesianDataFrame1D;

    fn mul(self, rhs: &CartesianDataFrame1D) -> Self::Output {
        let mut result = CartesianDataFrame1D::new_from(&rhs.underlying_mesh, rhs.n_comp, rhs.n_ghost);
        for (i, vals) in result.data.iter_mut().enumerate(){
            *vals = self * rhs.data[i];
        }
        result
    }
}

impl std::ops::Add<&CartesianDataFrame1D> for &CartesianDataFrame1D {
    type Output = CartesianDataFrame1D;

    fn add(self, rhs: &CartesianDataFrame1D) -> Self::Output {
        let mut sum = CartesianDataFrame1D::new_from(&rhs.underlying_mesh, rhs.n_comp, rhs.n_ghost);
        for (i, vals) in sum.data.iter_mut().enumerate() {
            *vals = self.data[i] + rhs.data[i];
        }
        sum
    }
}

impl std::ops::Mul<&CartesianDataFrame1D> for &CartesianDataFrame1D {
    type Output = CartesianDataFrame1D;

    fn mul(self, rhs: &CartesianDataFrame1D) -> Self::Output {
        let mut result = CartesianDataFrame1D::new_from(&rhs.underlying_mesh, rhs.n_comp, rhs.n_ghost);
        for (i, vals) in result.data.iter_mut().enumerate(){
            *vals = self.data[i] * rhs.data[i];
        }
        result
    }
}




//impl DataFrame for CartesianDataFrame {}


/// mutable iterator for `CartesianDataFrame`.
pub struct CartesianDataFrameGhostIter<'a> {
    df: &'a mut CartesianDataFrame1D,
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
    type Item = ((isize, usize), &'a mut Real);

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
    fn is_valid_cell () {
        let m1 = CartesianMesh1D::new(vec![0.0], vec![6.0], vec![3]);
        let df = CartesianDataFrame1D::new_from(&m1, 1, 1);

        assert_eq!(df.p_is_valid_cell(0), false);
        assert_eq!(df.p_is_valid_cell(1), true);
        assert_eq!(df.p_is_valid_cell(2), true);
        assert_eq!(df.p_is_valid_cell(3), true);
        assert_eq!(df.p_is_valid_cell(4), false);
    }

    #[test]
    fn enumerate_pos () {
        let m1 = CartesianMesh1D::new(vec![0.0], vec![6.0], vec![3]);
        let mut df = CartesianDataFrame1D::new_from(&m1, 1, 1);
        

        let mut df_iter_pos = df.into_iter().enumerate_pos();

        assert_eq!(df_iter_pos.next(), Some(((1.0, 0), &mut 0.0)));
        assert_eq!(df_iter_pos.next(), Some(((3.0, 0), &mut 0.0)));
        

    }

    #[test]
    // currently only tests the iterator on 1D data
    fn data_frame_iterator() {
        // 1D
        {
            let m1 = CartesianMesh1D::new(vec![0.0], vec![6.0], vec![3]);
            let mut df = CartesianDataFrame1D::new_from(&m1, 1, 1);
            df.fill_ic(|x,_|  x + 1.0);

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
        
    }


    #[test]
    fn mesh_iterator () {
        let m = CartesianMesh1D::new(vec![0.0], vec![10.0], vec![5]);
        let mut m_iter = m.into_iter();
        assert_eq!(m_iter.next().unwrap(), vec![&mut 1.0]);
        assert_eq!(m_iter.next().unwrap(), vec![&mut 3.0]);
        assert_eq!(m_iter.next().unwrap(), vec![&mut 5.0]);
        assert_eq!(m_iter.next().unwrap(), vec![&mut 7.0]);
        assert_eq!(m_iter.next().unwrap(), vec![&mut 9.0]);
        assert_eq!(m_iter.next(), Option::None);

    }

    #[test]
    fn indexing () {
        // 1D cartesian mesh
        let m1 = CartesianMesh1D::new(vec![0.0], vec![10.0], vec![5]);
        assert_eq!(m1.get_ijk(3).unwrap(), (3,0));
        assert_eq!(m1.get_ijk(5), Option::None);
        assert_eq!(m1.index(3, 0), &7.0);
        let mut cdf1 = CartesianDataFrame1D::new_from(&m1, 1, 2);
        cdf1.fill_ic(|x,_|x + 1.0);


    }

    #[test]
    fn node_pos () {
        let m1 = CartesianMesh1D::new(vec![0.0], vec![10.0], vec![5]);
        assert_eq!(m1.node_pos, vec![1.0, 3.0, 5.0, 7.0, 9.0]);
    }

    #[test]
    fn test_dirichlet_bc() {
        let u1 = CartesianMesh1D::new(vec![0.0], vec![10.0], vec![5]);
        let mut cdf = CartesianDataFrame1D::new_from(&u1, 1, 2);
        cdf.fill_ic(|x,_| x + 1.0);
        let bc = BCs::new(vec![
            ComponentBCs::new(
                    vec![BCType::Dirichlet(0.0)], 
                    vec![BCType::Dirichlet(1.0)],
            )]
        );
        cdf.fill_bc(&bc);
        assert_eq!(
            cdf.data,
            vec![-6.0, -2.0, 2.0, 4.0, 6.0, 8.0, 10.0, -8.0, -26.0]
        );
    }

    #[test]
    fn test_neumann_bc() {
        let u1 = CartesianMesh1D::new(vec![0.0], vec![10.0], vec![5]);
        let mut cdf = CartesianDataFrame1D::new_from(&u1, 1, 2);
        cdf.fill_ic(|x,_| x + 1.0);
        let bc = BCs::new(vec![
            ComponentBCs::new(
                    vec![BCType::Neumann(0.0)], 
                    vec![BCType::Neumann(1.0)],
            )]
        );
        cdf.fill_bc(&bc);
        assert_eq!(
            cdf.data,
            vec![4.0, 2.0, 2.0, 4.0, 6.0, 8.0, 10.0, 10.0, 8.0]
        );
    }

    #[test]
    fn test_prescribed_bc() {
        let u1 = CartesianMesh1D::new(vec![0.0], vec![10.0], vec![5]);
        let mut cdf = CartesianDataFrame1D::new_from(&u1, 1, 2);
        cdf.fill_ic(|x,_| x + 1.0);
        let bc = BCs::new(vec![
            ComponentBCs::new(
                    vec![BCType::Prescribed(vec![-1.0, -2.0])], 
                    vec![BCType::Prescribed(vec![15.0, 16.0])],
            )]
        );
        cdf.fill_bc(&bc);
        assert_eq!(
            cdf.data,
            vec![-2.0, -1.0, 2.0, 4.0, 6.0, 8.0, 10.0, 15.0, 16.0]
        );
    }
}

