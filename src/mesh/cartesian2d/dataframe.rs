
// ********************************* Cartesian data frame ****************************************
use std::rc::Rc;
//use crate::{Real, UIntVec2};
use super::*;
use super::mesh::*;
use crate::boundary_conditions::{BCs, BcType};
use super::super::DataFrame;

/// Enum to mark each cell as either valid, or the location of the ghost cell
#[derive(Debug, Clone, Eq, PartialEq)]
pub enum CellType{
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
pub struct CartesianDataFrame2D<T>{
    /// The data is stored here
    pub data: Vec<T>,

    /// The number of ghost nodes, added to each side of the underlying `CartesianMesh`
    pub n_ghost: usize,

    /// The number of cells in each direction, after the ghost cells have been added
    pub n_grown: UIntVec2,

    /// Reference to the underlying `CartesianMesh`
    pub underlying_mesh: Rc<CartesianMesh2D>,

    /// The number of components to be stored in the data frame. For example, storing
    /// a 3D velocity would result in `n_comp = 3`
    pub n_comp: usize,

    /// The total number of individual pieces of information needed to be stored
    pub n_nodes: usize,

    /// Distinguish between valid and ghost cells
    cell_type: Vec<CellType>,

    /// The boundary conditions to apply for the data frame
    pub bc: BCs,
}

impl <T: Clone + Default> DataFrame for CartesianDataFrame2D<T> {}

/// data structure to store data on CartesianMesh
impl <T: Clone+Default> CartesianDataFrame2D<T> {
    /// generate new `CartesianDataFrame` from a given mesh, adding a given
    /// number of ghost nodes
    pub fn new_from(m: & Rc<CartesianMesh2D>, 
                    bc: BCs,
                    n_comp: usize, 
                    n_ghost: usize) -> CartesianDataFrame2D<T>
    {
        let n_nodes = ((m.n[0] + 2*n_ghost) * (m.n[1] + 2*n_ghost))*n_comp;

        // calculate grown size of the data frame
        let n_grown: UIntVec2 = m.n.clone().iter().map(|n| n + 2*n_ghost).collect();

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
                if j < 0 { *cell = CellType::SouthWest; }
                else if j >= m.n[1] as isize { *cell = CellType::NorthWest; }
                else { *cell = CellType::West; }
            }
            else if i >= m.n[0] as isize {
                if j < 0 { *cell = CellType::SouthEast; }
                else if j >= m.n[0] as isize { *cell = CellType::NorthEast; } 
                else{ *cell = CellType::East; }
            }
        }

        CartesianDataFrame2D
        {
            n_ghost,
            data: vec![std::default::Default::default(); n_nodes],
            n_grown,
            underlying_mesh: Rc::clone(m),
            n_comp,
            n_nodes,
            cell_type,
            bc,
        }
    }

    /// Fill `CartesianDataFrame` from a initial condition function
    pub fn fill_ic (&mut self, ic: impl Fn(Real, Real, usize)->T)
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
    pub fn get_valid_cells(&self) -> Vec<T> {
        self.data
            .clone()
            .into_iter()
            .enumerate()
            .filter(|&(p, _)| self.p_is_valid_cell(p))
            .map(|(_,val)| val)
            .collect()
    }

    /// create an immutable iterator over the data
    pub fn iter(&self) -> CartesianDataFrame2DIter<T> {
        CartesianDataFrame2DIter{
            df: self,
            current_indx: 0,
        }
    } 

    /// create a mutable iterator over the data
    pub fn iter_mut(&mut self) -> CartesianDataFrame2DIterMut<T> {
        CartesianDataFrame2DIterMut {
            df: self,
            current_indx: 0,
        }
    }
}

impl <T> core::ops::IndexMut<(isize, isize, usize)> for CartesianDataFrame2D<T> {
    /// Exclusively borrows element at (i,j,n). The valid cells are indexed from zero 
    /// and ghost cells at the lower side of the domain are indexed with negative numbers
    fn index_mut(&mut self, indx: (isize,isize,usize)) -> &mut T {
        let (i,j,n) = indx;
        let p = get_flat_index_unchecked(i,j,n,&self.n_grown,self.n_ghost);
        &mut self.data[p]
    }
}


impl <T> core::ops::Index<(isize, isize, usize)> for CartesianDataFrame2D<T>{
    type Output = T;

    /// Borrows the element at (i,j,n). The valid cells are indexed from zero
    /// and ghost cells at the lower side of the domain are indexed with negative
    /// numbers.
    fn index(&self, indx: (isize, isize, usize) ) -> & Self::Output {
        let (i,j,n) = indx;
        let p = get_flat_index_unchecked(i, j, n, &self.n_grown, self.n_ghost);

        &self.data[p]
    }
}


pub trait BoundaryCondition2D<T> {
    /// function which fills the boundary conditions
    fn fill_bc(&mut self);

    /// gather ghost cell data into a vector to move to another block
    fn collect_ghost_cells(&self, side: Vec<CellType>) -> Vec<T>;
}

trait GhostCells{
    /// checks if the cell at (i,j,k) contains a valid or a ghost cell. Returns true if valid,
    /// and returns false if ghost
    fn ij_is_valid_cell(&self, i: isize, j: isize) -> bool;

    /// check if the cell at p contains a valid or ghost cell. Returns same as ijk_is_valid_cell
    fn p_is_valid_cell(&self, p: usize) -> bool;
}

// private implementation of boundary conditions
impl CartesianDataFrame2D<Real> {
    fn fill_neumann_low(&mut self, gradient: &Real, i_comp: usize, i_dim: usize) {
        for i in 0isize..self.underlying_mesh.n[i_dim] as isize{
            if i_dim == 0 {
                for ghost_indx in 1..self.n_ghost as isize + 1{
                    self[(-ghost_indx, i,i_comp)] = self[(ghost_indx-1,i,i_comp)] 
                            - gradient * (2.0*ghost_indx as Real - 1.0)*self.underlying_mesh.dx[i_dim];
                }
            }
            if i_dim == 1 {
                for ghost_indx in 1..self.n_ghost as isize + 1{
                    self[(i, -ghost_indx,i_comp)] = self[(i,ghost_indx-1,i_comp)] 
                            - gradient * (2.0*ghost_indx as Real-1.0)*self.underlying_mesh.dx[i_dim];
                }
            }
        }
    }

    fn fill_neumann_high(&mut self, gradient: &Real, i_comp: usize, i_dim: usize) {
        let hi = [self.underlying_mesh.n[0] as isize, self.underlying_mesh.n[1] as isize];
        for i in 0isize..self.underlying_mesh.n[i_dim] as isize{
            if i_dim == 0{
                for ghost_indx in 0..self.n_ghost as isize{
                    self[(hi[i_dim]+ghost_indx, i,i_comp)] = self[(hi[i_dim]-1-ghost_indx,i,i_comp)] 
                      + gradient * (2.0*(ghost_indx+1) as Real-1.0)*self.underlying_mesh.dx[i_dim];
                }
            }
            if i_dim == 1 {
                for ghost_indx in 0..self.n_ghost as isize {
                    self[(i,hi[i_dim]+ghost_indx,i_comp)] = self[(i,hi[i_dim]-1-ghost_indx,i_comp)] 
                      + gradient * (2.0*(ghost_indx+1) as Real-1.0)*self.underlying_mesh.dx[i_dim];
                }
            }
        }
    }

    fn fill_dirichlet_low(&mut self, value: &Real, i_comp: usize, i_dim: usize) {
        for i in 0isize..self.underlying_mesh.n[i_dim] as isize{
            if i_dim == 0{
                for ghost_indx in 0..self.n_ghost as isize{
                    self[(-ghost_indx-1,i,i_comp)] = 2.0*value - self[(ghost_indx,i,i_comp)];
                }
            }
            else if i_dim == 1 {
                for ghost_indx in 0..self.n_ghost as isize{
                    self[(i,-ghost_indx-1,i_comp)] = 2.0*value - self[(i,ghost_indx,i_comp)];
                }
            }
        }
    }

    fn fill_dirichlet_high(&mut self, value: &Real, i_comp: usize, i_dim: usize) {
        let hi = vec![self.underlying_mesh.n[0] as isize, self.underlying_mesh.n[1] as isize];
        for i in 0isize..self.underlying_mesh.n[i_dim] as isize{
            if i_dim == 0{
                for ghost_indx in 0..self.n_ghost as isize{
                    self[(hi[i_dim]+ghost_indx,i,i_comp)] = 2.0*value 
                                                - self[(hi[i_dim]-1-ghost_indx,i,i_comp)];
                }
            }
            else if i_dim == 1 {
                for ghost_indx in 0..self.n_ghost as isize{
                    self[(i,hi[i_dim]+ghost_indx,i_comp)] = 2.0*value 
                                                - self[(i,hi[i_dim]-1-ghost_indx,i_comp)];
                }
            }
        }
    }
}

impl <'a> BoundaryCondition2D<Real> for CartesianDataFrame2D<Real>{
    /// Fill the ghost nodes of the CartesianDataFrame based on BcType
    fn fill_bc (&mut self) {
        for i_comp in 0..self.n_comp{
            for i_dim in 0..2{
                // cloning these is a dirty hack to avoid issues with borrowing
                let bc_lo = &self.bc.bcs[i_comp].lo[i_dim].clone();
                let bc_hi = &self.bc.bcs[i_comp].hi[i_dim].clone();
                // low boundary condition
                match bc_lo {
                    BcType::Neumann(gradient) => {
                        self.fill_neumann_low(gradient, i_comp, i_dim);
                    }
                    BcType::Dirichlet(value) => {
                        self.fill_dirichlet_low(value, i_comp, i_dim);
                    }
                    BcType::Reflect => {
                        panic!("Reflect boundary condition not supported yet")
                    }
                    BcType::Internal(_id) => {
                        panic!("Internal boundary condition not supported yet");
                    }
                }

                // high boundary condition
                match bc_hi {
                    BcType::Neumann(gradient) => {
                        self.fill_neumann_high(gradient, i_comp, i_dim);
                    }
                    
                    BcType::Dirichlet(value) => {
                        self.fill_dirichlet_high(value, i_comp, i_dim);
                    }
                    BcType::Reflect => {
                        panic!("Reflect boundary condition not supported yet")
                    }
                    BcType::Internal(_id) => {
                        panic!("Internal boundary condition not supported yet")
                    }
                }
            }
        }
    }

    fn collect_ghost_cells(&self, side: Vec<CellType>) -> Vec<Real> {
        self.iter()
            .ghost(side)
            .cloned()
            .collect()
    }
}   

impl<T: Clone+Default> GhostCells for CartesianDataFrame2D<T>{
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

/// Immutable iterator for `CartesianDataFrame`.
pub struct CartesianDataFrame2DIter<'a, T> {
    df: &'a CartesianDataFrame2D<T>,
    current_indx: usize,
}

impl <'a, T> Iterator for CartesianDataFrame2DIter<'a, T>{
    type Item = &'a T;

    fn next(&mut self) -> Option<Self::Item> {
        // progress the current index to skip ghost cells
        while self.current_indx <= self.df.n_nodes - self.df.n_comp &&
                                  // !self.df.p_is_valid_cell(self.current_indx) {
                                  !(self.df.cell_type[self.current_indx] == CellType::Valid) {
            self.current_indx += 1;
        }

        // check if the next item exists
        if self.current_indx <= self.df.n_nodes - self.df.n_comp
                   //&& self.df.p_is_valid_cell(self.current_indx){
                    && self.df.cell_type[self.current_indx] == CellType::Valid {

            // access and return next item
            let next_data = Some(&self.df.data[self.current_indx]);
            self.current_indx += 1;
            next_data
        }
        else{
            // the next item doesn't exist, so return None
            None
        }
    }
}


/// Mutable iterator for `CartesianDataFrame`.
pub struct CartesianDataFrame2DIterMut<'a, T> {
    df: &'a mut CartesianDataFrame2D<T>,
    current_indx: usize,
}
impl <'a, T> Iterator for CartesianDataFrame2DIterMut<'a, T> {
    type Item = &'a mut T;

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


impl <'a, T> IntoIterator for &'a mut CartesianDataFrame2D<T>{
    type Item = &'a mut T;
    type IntoIter = CartesianDataFrame2DIterMut<'a,T>;

    fn into_iter(self) -> Self::IntoIter{
        Self::IntoIter {
            current_indx: 0,
            df: self,
        }
    }
}

impl<'a, T> IntoIterator for &'a CartesianDataFrame2D<T>{
    type Item = &'a T;
    type IntoIter = CartesianDataFrame2DIter<'a, T>;

    fn into_iter(self) -> Self::IntoIter{
        Self::IntoIter{
            current_indx: 0,
            df: self,
        }
    }
}


/// Struct to mutably iterate over a data frame, with the current index and the data itself
pub struct EnumerateIndex<'a, T>{
    iter: CartesianDataFrame2DIterMut<'a, T>,
}

/// Turns `CartesianDataFrameIter` into a mutable iterator which enumerates over the current index
pub trait IndexEnumerable <'a, T> {
    fn enumerate_index(self) -> EnumerateIndex<'a, T>;
}
impl <'a,T> IndexEnumerable <'a, T> for CartesianDataFrame2DIterMut<'a,T> {
    fn enumerate_index(self) -> EnumerateIndex<'a,T> {
        EnumerateIndex{
            iter: self,
        }
    }
}

impl <'a,T: Clone+Default> Iterator for EnumerateIndex<'a,T> {
    type Item = ((isize, isize, usize), &'a mut T);

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

/// Structure to mutably iterate over data, giving both the position of the data,
///  and the data itself
pub struct EnumeratePos<'a,T>{
    iter: CartesianDataFrame2DIterMut<'a,T>,
}

/// Turns `CartesianDataFrameIter` into a mutable iterator which enumerates the position 
/// of the data
pub trait PosEnumerable <'a,T> {
    fn enumerate_pos(self) -> EnumeratePos<'a,T>;
}
impl <'a,T> PosEnumerable <'a,T> for CartesianDataFrame2DIterMut<'a,T> {
    fn enumerate_pos(self) -> EnumeratePos<'a,T> {
        EnumeratePos{
            iter: self,
        }
    }
}

impl <'a,T: Clone+Default> Iterator for EnumeratePos<'a,T> {
    type Item = ((Real, Real, usize), &'a mut T);

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

// multiplication of owned data frame by a real number
impl std::ops::Mul<CartesianDataFrame2D<Real>> for Real {
    type Output = CartesianDataFrame2D<Real>;

    fn mul(self, rhs: CartesianDataFrame2D<Real>) -> Self::Output {
        let mut result = CartesianDataFrame2D::new_from(
            &rhs.underlying_mesh, rhs.bc.clone(), rhs.n_comp, rhs.n_ghost
        );
        for (i, vals) in result.data.iter_mut().enumerate(){
            *vals = self * rhs.data[i];
        }
        result
    }
}

/// multiplication of borrowed data frame by real
/// returns a new owned data frame
/// takes the boundary conditions of the rhs data frame
impl std::ops::Mul<&CartesianDataFrame2D<Real>> for Real {
    type Output = CartesianDataFrame2D<Real>;

    fn mul(self, rhs: &CartesianDataFrame2D<Real>) -> Self::Output {
        let mut result = CartesianDataFrame2D::new_from(
            &rhs.underlying_mesh, rhs.bc.clone(), rhs.n_comp, rhs.n_ghost
        );
        for (i, vals) in result.data.iter_mut().enumerate(){
            *vals = self * rhs.data[i];
        }
        result
    }
}

/// addition of two borrowed data frames
/// returns a new owned data frame
/// takes the boundary conditions of the rhs data frame
impl std::ops::Add<&CartesianDataFrame2D<Real>> for &CartesianDataFrame2D<Real> {
    type Output = CartesianDataFrame2D<Real>;

    fn add(self, rhs: &CartesianDataFrame2D<Real>) -> Self::Output {
        let mut sum = CartesianDataFrame2D::new_from(
            &rhs.underlying_mesh, rhs.bc.clone(), rhs.n_comp, rhs.n_ghost
        );
        for (i, vals) in sum.data.iter_mut().enumerate() {
            *vals = self.data[i] + rhs.data[i];
        }
        sum
    }
}

/// multiplication of two borrowed data frames
/// returns a new owned data frame
/// takes the boundary conditions of the rhs data frame
impl std::ops::Mul<&CartesianDataFrame2D<Real>> for &CartesianDataFrame2D<Real> {
    type Output = CartesianDataFrame2D<Real>;

    fn mul(self, rhs: &CartesianDataFrame2D<Real>) -> Self::Output {
        let mut result = CartesianDataFrame2D::new_from(
            &rhs.underlying_mesh, rhs.bc.clone(), rhs.n_comp, rhs.n_ghost
        );
        for (i, vals) in result.data.iter_mut().enumerate(){
            *vals = self.data[i] * rhs.data[i];
        }
        result
    }
}

pub struct CartesianDataFrameGhostIter<'a, T> {
    df: &'a CartesianDataFrame2D<T>,
    current_indx: usize,
    side: Vec<CellType>,
}

/// mutable iterator over the east cells in `CartesianDataFrame`.
pub struct CartesianDataFrameGhostIterMut<'a,T> {
    df: &'a mut CartesianDataFrame2D<T>,
    current_indx: usize,
    side: Vec<CellType>,
}

impl <'a,T: Clone + Default> Iterator for CartesianDataFrameGhostIterMut<'a,T> {
    type Item = &'a mut T;

    // this is a safe function wrapping unsafe code. Rust cannot guarantee that it is safe, but
    // in practice it can be, and I'm pretty sure that it should be safe for everything we will 
    // ever want to do
    fn next(&mut self) -> Option<Self::Item> {
        // progress the current index to next on the given side
        while self.current_indx <= self.df.n_nodes - self.df.n_comp &&
                                  !(self.side.contains(&self.df.cell_type[self.current_indx])) {
            self.current_indx += 1;
        }

        // check if the next item exists
        if self.current_indx <= self.df.n_nodes-self.df.n_comp {
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

impl <'a, T: Clone + Default> Iterator for CartesianDataFrameGhostIter<'a, T> {
    type Item = &'a T;

    fn next(&mut self) -> Option<Self::Item> {
        // progress the current index to next on the given side
        while self.current_indx <= self.df.n_nodes - self.df.n_comp &&
                                  !(self.side.contains(&self.df.cell_type[self.current_indx])) {
            self.current_indx += 1;
        }

        // check if the next item exists
        if self.current_indx <= self.df.n_nodes-self.df.n_comp{
            // access and return next item
            let next_data = Some(&self.df.data[self.current_indx]);
            self.current_indx += 1;
            next_data
        }
        // if the next item doesn't exist, return None
        else {
            None
        }
    }
}

pub trait GhostIteratorMut <'a,T> {
    fn ghost(self, side: Vec<CellType>) -> CartesianDataFrameGhostIterMut<'a,T>;
}

pub trait GhostIterator <'a, T> {
    fn ghost(self, side: Vec<CellType>) -> CartesianDataFrameGhostIter<'a, T>;
}

impl <'a,T> GhostIteratorMut <'a,T> for CartesianDataFrame2DIterMut<'a,T> {
    fn ghost(self, side: Vec<CellType>) -> CartesianDataFrameGhostIterMut<'a,T> {
        CartesianDataFrameGhostIterMut{
            current_indx: 0,
            df: self.df,
            side,
        }   
    }
}

impl <'a, T> GhostIterator <'a, T> for CartesianDataFrame2DIter<'a, T> {
    fn ghost(self, side: Vec<CellType>) -> CartesianDataFrameGhostIter<'a, T> {
        CartesianDataFrameGhostIter{
            current_indx: 0,
            df: self.df,
            side,
        }
    }
}


/// Struct to iterate over a data frame, with the current index and the data itself
pub struct EnumerateGhostIndex<'a,T>{
    iter: CartesianDataFrameGhostIterMut<'a,T>,
}

/// Turns `CartesianDataFrameIter` into an iterator which enumerates over the current index
pub trait GhostIndexEnumerable <'a,T> {
    fn enumerate_index(self) -> EnumerateGhostIndex<'a,T>;
}

impl <'a,T> GhostIndexEnumerable <'a,T> for CartesianDataFrameGhostIterMut<'a,T> {
    fn enumerate_index(self) -> EnumerateGhostIndex<'a,T> {
        EnumerateGhostIndex{
            iter: self,
        }
    }
}

impl <'a,T: Clone + Default> Iterator for EnumerateGhostIndex<'a,T> {
    type Item = ((isize, isize, usize), &'a mut T);

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
    use crate::boundary_conditions::*;

    #[test]
    fn data_frame_iterator() {
        
        // 2D
        {
            let m2 = CartesianMesh2D::new(
                RealVec2([0.0, 0.0]), RealVec2([6.0, 6.0]), UIntVec2([3, 3])
            );
            let mut df2 = CartesianDataFrame2D::<Real>::new_from(&m2, BCs::empty(), 1, 1);
            df2.fill_ic(|x,y,_| x + y);
            
            let mut df2_iter = df2.into_iter();
            assert_eq!(df2_iter.next(), Some(&2.0));
            assert_eq!(df2_iter.next(), Some(& 4.0));
            assert_eq!(df2_iter.next(), Some(& 6.0));
            assert_eq!(df2_iter.next(), Some(& 4.0));
            assert_eq!(df2_iter.next(), Some(& 6.0));
            assert_eq!(df2_iter.next(), Some(& 8.0));
            assert_eq!(df2_iter.next(), Some(& 6.0));
            assert_eq!(df2_iter.next(), Some(& 8.0));
            assert_eq!(df2_iter.next(), Some(& 10.0));
            assert_eq!(df2_iter.next(), None);
        }
    }

    #[test]
    fn node_type () {
        let m2 = CartesianMesh2D::new(RealVec2([0.0,0.0]), RealVec2([6.0, 6.0]), UIntVec2([3,3]));
        let df2 = CartesianDataFrame2D::<Real>::new_from(&m2, BCs::empty(), 1, 1);

        assert_eq!(
            df2.cell_type,
            vec![
                CellType::SouthWest,
                CellType::South,
                CellType::South,
                CellType::South,
                CellType::SouthEast,
                CellType::West,
                CellType::Valid,
                CellType::Valid,
                CellType::Valid,
                CellType::East,
                CellType::West,
                CellType::Valid,
                CellType::Valid,
                CellType::Valid,
                CellType::East,
                CellType::West,
                CellType::Valid,
                CellType::Valid,
                CellType::Valid,
                CellType::East,
                CellType::NorthWest,
                CellType::North,
                CellType::North,
                CellType::North,
                CellType::NorthEast
            ]
        );

    }

    #[test]
    fn test_dirichlet_bc_one_ghost() {
        let u1 = CartesianMesh2D::new(RealVec2([0.0, 0.0]), RealVec2([6.0, 6.0]), UIntVec2([3,3]));
        let bc = BCs::new(vec![
            ComponentBCs::new(
                //        x direction             y direction
                vec![BcType::Dirichlet(0.0), BcType::Dirichlet(3.0)], // low
                vec![BcType::Dirichlet(7.0), BcType::Dirichlet(4.0)], //high
            )
        ]);
        let mut cdf = CartesianDataFrame2D::new_from(&u1, bc, 1, 1);
        cdf.fill_ic(|x,_,_| x + 1.0);
        cdf.fill_bc();
        assert_eq!(
            cdf.data,
            vec![
                 0.0, 4.0, 2.0, 0.0, 0.0, 
                -2.0, 2.0, 4.0, 6.0, 8.0,
                -2.0, 2.0, 4.0, 6.0, 8.0,
                -2.0, 2.0, 4.0, 6.0, 8.0,
                 0.0, 6.0, 4.0, 2.0, 0.0     
            ]
        );
    }

    #[test]
    fn test_neumann_bc_one_ghost() {
        let u1 = CartesianMesh2D::new(RealVec2([0.0, 0.0]), RealVec2([6.0, 6.0]), UIntVec2([3,3]));
        let bc = BCs::new(vec![
            ComponentBCs::new(
                //        x direction             y direction
                vec![BcType::Neumann( 0.0), BcType::Neumann(-1.0)], // low
                vec![BcType::Neumann(-2.0), BcType::Neumann( 1.0)], // high
            )
        ]);
        let mut cdf = CartesianDataFrame2D::new_from(&u1, bc, 1, 1);
        cdf.fill_ic(|x,_,_| x + 1.0);
        cdf.fill_bc();
        assert_eq!(
            cdf.data,
            vec![
                0.0, 4.0, 6.0, 8.0, 0.0,
                2.0, 2.0, 4.0, 6.0, 2.0, 
                2.0, 2.0, 4.0, 6.0, 2.0,
                2.0, 2.0, 4.0, 6.0, 2.0,
                0.0, 4.0, 6.0, 8.0, 0.0
            ]
        )
    }

    #[test]
    fn test_neumann_bc_two_ghost(){
        let u1 = CartesianMesh2D::new(RealVec2([0.0, 0.0]), RealVec2([6.0, 6.0]), UIntVec2([3,3]));
        let bc = BCs::new(vec![
            ComponentBCs::new(
                //        x direction             y direction
                vec![BcType::Neumann( 0.0), BcType::Neumann(-1.0)], // low
                vec![BcType::Neumann(-2.0), BcType::Neumann( 1.0)], // high
            )
        ]);
        let mut cdf = CartesianDataFrame2D::new_from(&u1, bc, 1, 2);
        cdf.fill_ic(|x,_,_| x + 1.0);
        cdf.fill_bc();
        assert_eq!(
            cdf.data,
            vec![
                0.0, 0.0, 8.0, 10.0, 12.0, 0.0,  0.0, 
                0.0, 0.0, 4.0,  6.0,  8.0, 0.0,  0.0,
                4.0, 2.0, 2.0,  4.0,  6.0, 2.0, -8.0,
                4.0, 2.0, 2.0,  4.0,  6.0, 2.0, -8.0,
                4.0, 2.0, 2.0,  4.0,  6.0, 2.0, -8.0,
                0.0, 0.0, 4.0,  6.0,  8.0, 0.0,  0.0,
                0.0, 0.0, 8.0, 10.0, 12.0, 0.0,  0.0
            ]
        )
    }

    #[test]
    fn test_dirichlet_bc_two_ghost () {
        let u1 = CartesianMesh2D::new(RealVec2([0.0, 0.0]), RealVec2([6.0, 6.0]), UIntVec2([3,3]));
        let bc = BCs::new(vec![
            ComponentBCs::new(
                //        x direction             y direction
                vec![BcType::Dirichlet(0.0), BcType::Dirichlet(3.0)], // low
                vec![BcType::Dirichlet(7.0), BcType::Dirichlet(4.0)], //high
            )
        ]);
        let mut cdf = CartesianDataFrame2D::new_from(&u1, bc, 1, 2);
        cdf.fill_ic(|x,_,_| x + 1.0);
        cdf.fill_bc();
        assert_eq!(
            cdf.data,
            vec![
                 0.0,  0.0, 4.0, 2.0, 0.0, 0.0,  0.0,
                 0.0,  0.0, 4.0, 2.0, 0.0, 0.0,  0.0,
                -4.0, -2.0, 2.0, 4.0, 6.0, 8.0, 10.0,
                -4.0, -2.0, 2.0, 4.0, 6.0, 8.0, 10.0,
                -4.0, -2.0, 2.0, 4.0, 6.0, 8.0, 10.0,
                 0.0,  0.0, 6.0, 4.0, 2.0, 0.0,  0.0,
                 0.0,  0.0, 6.0, 4.0, 2.0, 0.0,  0.0
            ]
        )
    }

    #[test]
    fn test_ghost_iter_two_ghost () {
        let u1 = CartesianMesh2D::new(RealVec2([0.0, 0.0]), RealVec2([6.0, 6.0]), UIntVec2([3,3]));
        let bc = BCs::new(vec![
            ComponentBCs::new(
                //        x direction             y direction
                vec![BcType::Dirichlet(0.0), BcType::Dirichlet(3.0)], // low
                vec![BcType::Dirichlet(7.0), BcType::Dirichlet(4.0)], //high
            )
        ]);
        let mut cdf = CartesianDataFrame2D::new_from(&u1, bc, 1, 2);
        cdf.fill_ic(|x,_,_| x + 1.0);
        cdf.fill_bc();
        //vec![
        //     0.0,  0.0, 4.0, 2.0, 0.0, 0.0,  0.0,
        //     0.0,  0.0, 4.0, 2.0, 0.0, 0.0,  0.0,
        //    -4.0, -2.0, 2.0, 4.0, 6.0, 8.0, 10.0,
        //    -4.0, -2.0, 2.0, 4.0, 6.0, 8.0, 10.0,
        //    -4.0, -2.0, 2.0, 4.0, 6.0, 8.0, 10.0,
        //     0.0,  0.0, 6.0, 4.0, 2.0, 0.0,  0.0,
        //     0.0,  0.0, 6.0, 4.0, 2.0, 0.0,  0.0
        //]
        let mut east_iter = cdf.iter_mut().ghost(vec![CellType::East]);
        assert_eq!(east_iter.next(), Some(&mut 8.0));
        assert_eq!(east_iter.next(), Some(&mut 10.0));
        assert_eq!(east_iter.next(), Some(&mut 8.0));
        assert_eq!(east_iter.next(), Some(&mut 10.0));
        assert_eq!(east_iter.next(), Some(&mut 8.0));
        assert_eq!(east_iter.next(), Some(&mut 10.0));
        assert_eq!(east_iter.next(), None);
    }
}
