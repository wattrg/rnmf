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
pub struct CartesianMesh2D {
    /// The position of the nodes in the mesh
    pub node_pos: Vec<Real>,

    /// The lower corner of the mesh
    pub lo: [Real; 2],

    /// The upper corner of the mesh
    pub hi: [Real; 2],

    /// The number of cells in nodes in each direction
    pub n : [usize; 2],

    /// The distance between each node in each
    pub dx: [Real; 2],

    pub dim: usize,
    
    /// The number of nodes in the mesh
    n_nodes: usize
}

/// provides functionality to assist with indexing of data. This is separate to the index operator
trait Indexing{
    /// Returns the non-flat index of the item stored at a particular flattened location
    fn get_ij(&self, p: usize ) -> Option<(isize, isize, usize)>;

    /// Returns a borrow to the value stored at a particular non-flat index
    fn index(&self, i: isize, j: isize, n: usize) -> &Real;
}


impl CartesianMesh2D {
    /// Generate new CartesianMesh from the lo and high corner, 
    /// and the number of nodes
    pub fn new(lo: [Real; 2], hi: [Real; 2], n: [usize; 2]) -> Rc<CartesianMesh2D>
    {
        let mut dx = [0.0; 2];
        for i_dim in 0..2 as usize{
            dx[i_dim] = (hi[i_dim] - lo[i_dim])/(n[i_dim] as Real);
        }

        // calculate the capacity of the vector required to store all the node positions
        let n_nodes = n[0]*n[1]*2;

        // allocate memory to the node positions vector
        let node_pos = Vec::with_capacity(n_nodes);

        // allocate memory to the mesh
        let mut cm = CartesianMesh2D
        {
            lo,
            hi,
            n,
            dx,
            node_pos,
            dim: 2,
            n_nodes,
        };

        // calculate the positions of the nodes
        // this will be able to be done better by creating an enumerating iterator
        // which will give (i,j,n) to begin with
        cm.node_pos = (0..n_nodes).map(
            |p: usize| -> Real {
                let (i,j,n) = cm.get_ij(p).unwrap();
                match n {
                    0 => {cm.lo[n] + ((i as Real) + 0.5) * cm.dx[n]}
                    1 => {cm.lo[n] + ((j as Real) + 0.5) * cm.dx[n]}
                    _ => {panic!("n cannot be {} in 2D", n)}
                }
                
            }).collect();
        Rc::new(cm)
    }
}

/// function to convert a multi-dimensional coordinate into a single dimension coordinate
fn get_flat_index(i: isize, j: isize, n: usize, 
            dimension: &[usize], ng: usize) -> usize {

    let mut p = i + ng as isize;
    p += dimension[0] as isize * (j + ng as isize) 
            + (dimension[0] * dimension[1] * n) as isize;
    

    p as usize
}

/// Converts (3+1)D index into a 1D (flattened) index
fn get_ij_from_p(p: usize, dimension: &[usize], n_nodes: usize, ng: usize) 
                                                        -> Option<(isize, isize, usize)>{
    if p >= n_nodes {
        Option::None
    }

    else{
        let i = (p % dimension[0]) as isize;
        let j;
        let d: isize;

        j = (p as isize - i)/dimension[0] as isize % dimension[1] as isize;

        d = (p as isize - j * dimension[1] as isize - i)
                    /dimension[0] as isize % 2;
        Some((i - ng as isize, j - ng as isize, d as usize))
    }
}

impl Indexing for CartesianMesh2D
{
    /// Retrieves the element at (i,j,n)
    fn index(&self, i: isize, j: isize, n: usize) -> & Real
    {
        let p = get_flat_index(i, j, n, &self.n, 0);
        let np = self.node_pos.get(p as usize);
        match np{
            Some(np) => {np},
            None => panic!("Index ({}, {}, {}) out of range of Cartesian mesh with size {:?}", 
                                i,j,n, self.n),
        }

    }

    /// Returns the un-flattened index
    fn get_ij(& self, p: usize) -> Option<(isize, isize, usize)>{
        get_ij_from_p(p, &self.n, self.n_nodes, 0)
    }
}


/// Immutable Iterator struct for CartesianMesh
pub struct CartesianMeshIter <'a> {
    mesh: &'a CartesianMesh2D,
    current_indx: usize,
}

impl <'a> Iterator for CartesianMeshIter<'a>{
    // return a vector of references instead of slice because the length of the returned object
    // depends on the number of dimensions, and thus is not known at compile time, so slices
    // won't work
    type Item = Vec<&'a Real>;
    
    fn next(& mut self) -> Option<Self::Item>{
        // If the next position doesn't exist, return None
        if 2*self.current_indx >= self.mesh.node_pos.len(){
            None
        }
        // If the next position exists, return a vector of size 2 containing references
        // to the position
        else {
            let mut next_pos = Vec::with_capacity(self.mesh.dim);
            let (i,j,n) = self.mesh.get_ij(self.current_indx).unwrap();
            for i_dim in 0..2 {
                next_pos.push(self.mesh.index(i, j, n + i_dim))
            }
            self.current_indx += 1;
            Some(next_pos)
        }
    }
}

impl <'a> IntoIterator for &'a CartesianMesh2D {
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
#[derive(Debug, Clone, Eq, PartialEq)]
enum CellType{
    Valid,
    North,
    South,
    East,
    West,
    NorthEast,
    SouthEast,
    NorthWest,
    SouthWest,
}

/// Structure to store data defined on a `CartesianMesh`
#[derive(Debug, Clone)]
pub struct CartesianDataFrame2D{
    /// The data is stored here
    pub data: Vec<Real>,

    /// The number of ghost nodes, added to each side of the underlying `CartesianMesh`
    pub n_ghost: usize,

    pub n_grown: Vec<usize>,

    /// Reference to the underlying `CartesianMesh`
    pub underlying_mesh: Rc<CartesianMesh2D>,

    /// The number of components to be stored in the data frame. For example, storing
    /// a 3D velocity would result in `n_comp = 3`
    pub n_comp: usize,

    /// The total number of individual pieces of information needed to be stored
    pub n_nodes: usize,

    /// Distinguish between valid and ghost cells
    cell_type: Vec<CellType>,
}


/// data structure to store data on CartesianMesh
impl CartesianDataFrame2D {
    /// generate new `CartesianDataFrame` from a given mesh, adding a given
    /// number of ghost nodes
    pub fn new_from(m: & Rc<CartesianMesh2D>, 
                    n_comp: usize, n_ghost: usize) -> CartesianDataFrame2D
    {
        let n_nodes = ((m.n[0] + 2*n_ghost) * (m.n[1] + 2*n_ghost))*n_comp;

        // calculate grown size of the data frame
        let n_grown: Vec<usize> = m.n.clone().iter().map(|n| n + 2*n_ghost).collect();

        // calculate cell_type
        let mut cell_type = vec![CellType::Valid; n_nodes];
        for (p,cell) in cell_type.iter_mut().enumerate(){
            // get the index in index space
            let (i,j,_) = get_ij_from_p(p, &n_grown, n_nodes, n_ghost).unwrap();

            // check if the cell is a ghost cell, and fill in the appropriate type of ghost cell if
            // it is
            if j < 0 { *cell = CellType::South; }
            else if j >= m.n[1] as isize { *cell = CellType::North; }

            if i < 0 {
                if j < 0 { *cell = CellType::SouthEast; }
                else if j >= m.n[1] as isize { *cell = CellType::NorthEast; }
                else { *cell = CellType::East; }
            }
            else if i >= m.n[0] as isize {
                if j < 0 { *cell = CellType::SouthWest; }
                else if j >= m.n[0] as isize { *cell = CellType::NorthWest; } 
                else{ *cell = CellType::West; }
            }
        }

        CartesianDataFrame2D
        {
            n_ghost,
            data:  vec![0.0; n_nodes],
            n_grown,
            underlying_mesh: Rc::clone(m),
            n_comp,
            n_nodes,
            cell_type,
        }
    }

    /// Fill `CartesianDataFrame` from a initial condition function
    pub fn fill_ic (&mut self, ic: impl Fn(Real, Real, usize)->Real)
    {
        for (pos, val) in self.into_iter().enumerate_pos(){
            let (x,y,n) = pos;
            *val = ic(x,y,n);
        }

    }

    /// Convert the flat index of a data point into the (3+1)D coordinates
    fn get_ij(&self, p: usize) -> Option<(isize, isize, usize)> {
        get_ij_from_p(p, &self.n_grown, self.n_nodes, self.n_ghost)
    }

    /// returns a flat vector containing only the valid cells
    pub fn get_valid_cells(&self) -> Vec<Real> {
        self.data
            .clone()
            .into_iter()
            .enumerate()
            .filter(|&(p, _)| self.p_is_valid_cell(p))
            .map(|(_,val)| val)
            .collect()
    }
}

impl core::ops::IndexMut<(isize, isize, usize)> for CartesianDataFrame2D {
    /// Exclusively borrows element at (i,j,n). The valid cells are indexed from zero 
    /// and ghost cells at the lower side of the domain are indexed with negative numbers
    fn index_mut(&mut self, indx: (isize,isize,usize)) -> &mut Real {
        let (i,j,n) = indx;
        let p = get_flat_index(i,j,n,&self.n_grown,self.n_ghost);
        &mut self.data[p]
    }
}


impl core::ops::Index<(isize, isize, usize)> for CartesianDataFrame2D{
    type Output = Real;

    #[allow(dead_code)] 
    /// Borrows the element at (i,j,n). The valid cells are indexed from zero
    /// and ghost cells at the lower side of the domain are indexed with negative
    /// numbers.
    fn index(&self, indx: (isize, isize, usize) ) -> &Real {
        let (i,j,n) = indx;
        let p = get_flat_index(i, j, n, &self.n_grown, self.n_ghost);

        &self.data[p]
    }
}


pub trait BoundaryCondition2D {
    /// function which fills the boundary conditions
    fn fill_bc(&mut self, bcs: &BCs);

    /// checks if the cell at (i,j,k) contains a valid or a ghost cell. Returns true if valid,
    /// and returns false if ghost
    fn ij_is_valid_cell(&self, i: isize, j: isize) -> bool;

    /// check if the cell at p contains a valid or ghost cell. Returns same as ijk_is_valid_cell
    fn p_is_valid_cell(&self, p: usize) -> bool;
}


impl <'a> BoundaryCondition2D for CartesianDataFrame2D{
    /// Fill the ghost nodes of the CartesianDataFrame based on BCType
    fn fill_bc (&mut self, bc: &BCs) {
        for i_comp in 0..self.n_comp{
            for i_dim in 0..2{
                let bc_lo = &bc.bcs[i_comp].lo[i_dim];
                let bc_hi = &bc.bcs[i_comp].hi[i_dim];
                // low boundary condition
                match bc_lo {
                    BCType::Prescribed(_values) => {
                        panic!("Prescribed boundary conditions not supported yet");
                    }
                    BCType::Neumann(gradient) => {
                        if self.n_ghost >= 2 {panic!{"more than two ghost cells not supported for Neumann BC yet"};}
                        for i in 0isize..self.underlying_mesh.n[i_dim] as isize{
                            if i_dim == 0 {
                                self[(-1, i,i_comp)] = self[(1,i,i_comp)] - 2.0 * gradient * self.underlying_mesh.dx[i_dim];
                            }
                            if i_dim == 1 {
                                self[(i, -1,i_comp)] = self[(i,1,i_comp)] - 2.0 * gradient * self.underlying_mesh.dx[i_dim];
                            }
                        }
                    }
                    BCType::Dirichlet(value) => {
                        for i in 0isize..self.underlying_mesh.n[i_dim] as isize{
                            if i_dim == 0{
                                let m = -(value - self[(0,i,i_comp)]);
                                self[(-1,i,i_comp)] = self[(0,i,i_comp)] - 2.0 * m;
                            }
                            else if i_dim == 1 {
                                let m = -(value - self[(i,0,i_comp)]);
                                self[(i,-1,i_comp)] = self[(i,0,i_comp)] - 2.0 * m;
                            }
                        }
                    }
                }

                // high boundary condition
                match bc_hi {
                    BCType::Prescribed(_values) => {
                        panic!("Prescribed boundary condition not yet supported");
                    }
                
                    BCType::Neumann(gradient) => {
                        if self.n_ghost >= 2 {panic!{"more than two ghost cells not supported for Neumann BC yet"};}
                        let hi = vec![self.underlying_mesh.n[0] as isize, self.underlying_mesh.n[1] as isize];
                        for i in 0isize..self.underlying_mesh.n[i_dim] as isize{
                            if i_dim == 0{
                                self[(hi[i_dim], i,i_comp)] = self[(hi[i_dim]-2,i,i_comp)] + 2.0 * gradient * self.underlying_mesh.dx[i_dim];
                            }
                            if i_dim == 1 {
                                self[(i,hi[i_dim],i_comp)] = self[(i,hi[i_dim]-2,i_comp)] + 2.0 * gradient * self.underlying_mesh.dx[i_dim];
                            }
                        }
                    }
                    
                    BCType::Dirichlet(value) => {
                        let hi = vec![self.underlying_mesh.n[0] as isize, self.underlying_mesh.n[1] as isize];
                        for i in 0isize..self.underlying_mesh.n[i_dim] as isize{
                            if i_dim == 0{
                                let m = value - self[(hi[i_dim]-1,i,i_comp)];
                                self[(hi[i_dim],i,i_comp)] = self[(hi[i_dim]-1,i,i_comp)] + 2.0 * m;
                            }
                            else if i_dim == 1 {
                                let m = value - self[(i, hi[i_dim]-1,i_comp)];
                                self[(i,hi[i_dim],i_comp)] = self[(i,hi[i_dim]-1,i_comp)] + 2.0 * m;
                            }
                        }
                    }
                }
            }
        }
    }
    

    /// Determines if the cell at (i,j) is a valid cell (returns true) or a ghost
    /// cell (returns false)
    /// This will eventually be depricated in favour of the cell_type element in the data structure
    fn ij_is_valid_cell(&self, i: isize, j: isize) -> bool {
        let above_zero = i >= 0 && j >= 0;

        let mut below_max = i < self.underlying_mesh.n[0] as isize;
        below_max = below_max && j < self.underlying_mesh.n[1]  as isize;
        above_zero && below_max
    }

    fn p_is_valid_cell(&self, p: usize) -> bool {
        let (i,j,_) = self.get_ij(p).unwrap();
        self.ij_is_valid_cell(i, j)
    }
}



/// mutable iterator for `CartesianDataFrame`.
pub struct CartesianDataFrameIter<'a> {
    df: &'a mut CartesianDataFrame2D,
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
                                  // !self.df.p_is_valid_cell(self.current_indx) {
                                  !(self.df.cell_type[self.current_indx] == CellType::Valid) {
            self.current_indx += 1;
        }

        // check if the next item exists
        if self.current_indx <= self.df.n_nodes-self.df.n_comp 
                   //&& self.df.p_is_valid_cell(self.current_indx){
                    && self.df.cell_type[self.current_indx] == CellType::Valid {
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

// impl <'a> Iterator for &'a CartesianDataFrameIter<'a>{
//     type Item = &'a Real;

//     fn next(self) -> Option<Self::Item>{
//         while self.current_indx <= self.df.n_nodes - self.df.n_comp &&
//                                    !self.df.p_is_valid_cell(self.current_indx){
//             self.current_indx += 1;
//         }

//         let next_data: Option::<Self::Item>;
//         if self.current_indx <= self.df.n_nodes - self.df.n_comp && self.df.p_is_valid_cell(self.current_indx){
//             next_data = Some(&self.df.data[self.current_indx]);
//             self.current_indx += 1;
//             next_data
//         }
//         else{
//             None
//         }
//     }
// }

impl <'a> IntoIterator for &'a mut CartesianDataFrame2D{
    type Item = &'a mut Real;
    type IntoIter = CartesianDataFrameIter<'a>;

    fn into_iter(self) -> Self::IntoIter{
        Self::IntoIter {
            current_indx: 0,
            df: self,
        }
    }
}

// impl <'a> IntoIterator for &'a CartesianDataFrame2D{
//     type Item = &'a Real;
//     type IntoIter = CartesianDataFrameIter<'a>;

//     fn into_iter(self) -> Self::IntoIter{
//         Self::IntoIter{
//             current_indx: 0,
//             df: self,
//         }
//     }
// }

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
    type Item = ((isize, isize, usize), &'a mut Real);

    fn next(&mut self) -> Option<Self::Item>{
        let next_item = self.iter.next();
        match next_item{
            Some(nv) => {
                Some((self.iter.df.get_ij(self.iter.current_indx-1).unwrap(), nv))
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
    type Item = ((Real, Real, usize), &'a mut Real);

    fn next(&mut self) -> Option<Self::Item>{
        let next_item = self.iter.next();
        match next_item{
            Some(nv) => {
                let (i,j,n) = self.iter.df.get_ij(self.iter.current_indx-1).unwrap();
                let pos_x = *self.iter.df.underlying_mesh.index(i, j, 0);
                let pos_y;
                pos_y = *self.iter.df.underlying_mesh.index(i, j, 1);
                Some(((pos_x, pos_y, n), nv))
            }
            None => {
                None
            }
        }
        
    }
}

impl std::ops::Mul<CartesianDataFrame2D> for Real {
    type Output = CartesianDataFrame2D;

    fn mul(self, rhs: CartesianDataFrame2D) -> Self::Output {
        let mut result = CartesianDataFrame2D::new_from(&rhs.underlying_mesh, rhs.n_comp, rhs.n_ghost);
        for (i, vals) in result.data.iter_mut().enumerate(){
            *vals = self * rhs.data[i];
        }
        result
    }
}
impl std::ops::Mul<&CartesianDataFrame2D> for Real {
    type Output = CartesianDataFrame2D;

    fn mul(self, rhs: &CartesianDataFrame2D) -> Self::Output {
        let mut result = CartesianDataFrame2D::new_from(&rhs.underlying_mesh, rhs.n_comp, rhs.n_ghost);
        for (i, vals) in result.data.iter_mut().enumerate(){
            *vals = self * rhs.data[i];
        }
        result
    }
}

impl std::ops::Add<&CartesianDataFrame2D> for &CartesianDataFrame2D {
    type Output = CartesianDataFrame2D;

    fn add(self, rhs: &CartesianDataFrame2D) -> Self::Output {
        let mut sum = CartesianDataFrame2D::new_from(&rhs.underlying_mesh, rhs.n_comp, rhs.n_ghost);
        for (i, vals) in sum.data.iter_mut().enumerate() {
            *vals = self.data[i] + rhs.data[i];
        }
        sum
    }
}

impl std::ops::Mul<&CartesianDataFrame2D> for &CartesianDataFrame2D {
    type Output = CartesianDataFrame2D;

    fn mul(self, rhs: &CartesianDataFrame2D) -> Self::Output {
        let mut result = CartesianDataFrame2D::new_from(&rhs.underlying_mesh, rhs.n_comp, rhs.n_ghost);
        for (i, vals) in result.data.iter_mut().enumerate(){
            *vals = self.data[i] * rhs.data[i];
        }
        result
    }
}




//impl DataFrame for CartesianDataFrame {}


/// mutable iterator for `CartesianDataFrame`.
pub struct CartesianDataFrameGhostIter<'a> {
    df: &'a mut CartesianDataFrame2D,
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
    type Item = ((isize, isize, usize), &'a mut Real);

    fn next(&mut self) -> Option<Self::Item>{
        let next_item = self.iter.next();
        match next_item{
            Some(nv) => {
                Some((self.iter.df.get_ij(self.iter.current_indx-1).unwrap(), nv))
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
    // currently only tests the iterator on 1D data
    fn data_frame_iterator() {
        
        // 2D
        {
            let m2 = CartesianMesh2D::new([0.0, 0.0], [6.0, 6.0], [3, 3]);
            let mut df2 = CartesianDataFrame2D::new_from(&m2, 1, 1);
            df2.fill_ic(|x,y,_| x + y);
            


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
    fn ghost_iterator() {
        let m2 = CartesianMesh2D::new([0.0, 0.0], [10.0, 10.0], [5, 5]);
        let mut df = CartesianDataFrame2D::new_from(&m2, 1, 1);
        let df_iter = df.into_iter().ghost().enumerate_index();
        let mut count = 0;
        for (_, _) in df_iter{
            count += 1;
        }
        assert_eq!(count, 24);
    }

    #[test]
    fn mesh_iterator () {

        let m2 = CartesianMesh2D::new([0.0, 0.0], [6.0, 6.0], [3, 3]);
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
        assert_eq!(m2_iter.next(), Option::None);
    }

    #[test]
    fn indexing () {

        let m2 = CartesianMesh2D::new([0.0,0.0], [6.0,6.0], [3,3]);
        assert_eq!(m2.get_ij(16).unwrap(), (1,2,1));
        assert_eq!(m2.get_ij(40), Option::None);
        assert_eq!(m2.index(2,1,1), &3.0);

    }

    #[test]
    fn node_pos () {

        let m2 = CartesianMesh2D::new([0.0, 0.0], [10.0, 10.0], [5,5]);
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

    }

    #[test]
    fn node_type () {
        let m2 = CartesianMesh2D::new([0.0,0.0], [6.0, 6.0], [3,3]);
        let df2 = CartesianDataFrame2D::new_from(&m2, 1, 1);

        assert_eq!(
            df2.cell_type,
            vec![
                CellType::SouthEast,
                CellType::South,
                CellType::South,
                CellType::South,
                CellType::SouthWest,
                CellType::East,
                CellType::Valid,
                CellType::Valid,
                CellType::Valid,
                CellType::West,
                CellType::East,
                CellType::Valid,
                CellType::Valid,
                CellType::Valid,
                CellType::West,
                CellType::East,
                CellType::Valid,
                CellType::Valid,
                CellType::Valid,
                CellType::West,
                CellType::NorthEast,
                CellType::North,
                CellType::North,
                CellType::North,
                CellType::NorthWest
            ]
        );

    }


}

