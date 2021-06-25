
// ********************************* Cartesian data frame ****************************************
use std::rc::Rc;
//use crate::{Real, UIntVec2};
use super::*;
use super::mesh::*;
use crate::boundary_conditions::{BCs, BcType, BoundaryCondition};
use super::super::{DataFrame};
use log::warn;

#[derive(Debug, Clone, Eq, PartialEq)]
pub enum Loc {
    Regular,
    North,
    South,
    East,
    West,
    NorthEast,
    NorthWest,
    SouthWest,
    SouthEast,
}

/// # CellType
///
/// Marks each cell with its location.
#[derive(Debug, Clone, Eq, PartialEq)]
pub enum CellType{
    Valid(Loc),
    Ghost(Loc),
}

/// Structure to store data defined on a `CartesianMesh`
#[derive(Debug, Clone)]
pub struct CartesianDataFrame2D<S>{
    /// The data is stored here
    pub data: Vec<S>,

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
    pub bc: BCs<S>,
}

// Mark CartesianDataFrame as DataFrame
impl DataFrame for CartesianDataFrame2D<Real>{}

/// data structure to store data on CartesianMesh
impl <S> CartesianDataFrame2D<S>
    where S: Clone + Default
{
    /// generate new `CartesianDataFrame` from a given mesh, adding a given
    /// number of ghost nodes
    pub fn new_from(m: & Rc<CartesianMesh2D>, 
                    bc: BCs<S>,
                    n_comp: usize, 
                    n_ghost: usize) -> CartesianDataFrame2D<S>
    {
        //assert!(m.n[0] >= 2 * n_ghost && m.n[1] >= 2 * n_ghost);
        if !(m.n[0] >= 2 * n_ghost) || !(m.n[1] >= 2 * n_ghost) {
            warn!("Dataframe too thin. Make sure you have at least twice the number
                     of ghost cells in each direction");
        }

        let n_nodes = ((m.n[0] + 2*n_ghost) * (m.n[1] + 2*n_ghost))*n_comp;

        // calculate grown size of the data frame
        let n_grown: UIntVec2 = m.n.clone().iter().map(|n| n + 2*n_ghost).collect();

        // calculate cell_type
        let mut cell_type: Vec<CellType> = vec![CellType::Valid(Loc::Regular); n_nodes];
        for (p,cell) in cell_type.iter_mut().enumerate(){
            // get the index in index space
            let (i,j,_) = get_ij_from_p(p, &n_grown, n_nodes, n_ghost).unwrap();

            // check for ghost cells in corners
            if j < 0 && i < 0{
                *cell = CellType::Ghost(Loc::SouthWest);
            }
            else if j < 0 && i >= m.n[0] as isize {
                *cell = CellType::Ghost(Loc::SouthEast);
            }
            else if i < 0 && j >= m.n[1] as isize {
                *cell = CellType::Ghost(Loc::NorthWest);
            }
            else if i >= m.n[0] as isize && j >= m.n[1] as isize {
                *cell = CellType::Ghost(Loc::NorthEast);
            }

            // check for edge ghost cells
            else if i < 0 {
                *cell = CellType::Ghost(Loc::West);
            }
            else if j < 0 {
                *cell = CellType::Ghost(Loc::South);
            }
            else if i >= m.n[0] as isize {
                *cell = CellType::Ghost(Loc::East);
            }
            else if j >= m.n[1] as isize {
                *cell = CellType::Ghost(Loc::North);
            }

            // check for corner valid edge cells
            else if i < n_ghost as isize && j < n_ghost as isize {
                *cell = CellType::Valid(Loc::SouthWest);
            }
            else if i >= (m.n[0] - n_ghost) as isize && j < n_ghost as isize {
                *cell = CellType::Valid(Loc::SouthEast);
            }
            else if j >= (m.n[1] - n_ghost) as isize && i < n_ghost as isize {
                *cell = CellType::Valid(Loc::NorthWest);
            }
            else if j >= (m.n[1] - n_ghost) as isize && i >= (m.n[0] - n_ghost) as isize{
                *cell = CellType::Valid(Loc::NorthEast);
            }

            // finally check for valid edge cells
            else if i < n_ghost as isize{
                *cell = CellType::Valid(Loc::West);
            }
            else if i >= (m.n[0] - n_ghost) as isize{
                *cell = CellType::Valid(Loc::East);
            }
            else if j < n_ghost as isize {
                *cell = CellType::Valid(Loc::South);
            }
            else if j >= (m.n[1] - n_ghost) as isize {
                *cell = CellType::Valid(Loc::North);
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
    pub fn fill_ic (&mut self, ic: impl Fn(Real, Real, usize)->S)
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
    pub fn get_valid_cells(&self) -> Vec<S> {
        self.data
            .clone()
            .into_iter()
            .enumerate()
            .filter(|&(p, _)| self.p_is_valid_cell(p))
            .map(|(_,val)| val)
            .collect()
    }

    /// create an immutable iterator over the data
    pub fn iter(&self) -> CartesianDataFrame2DIter<S> {
        CartesianDataFrame2DIter{
            df: self,
            current_indx: 0,
        }
    } 

    /// create a mutable iterator over the data
    pub fn iter_mut(&mut self) -> CartesianDataFrame2DIterMut<S> {
        CartesianDataFrame2DIterMut {
            df: self,
            current_indx: 0,
        }
    }
}

impl <S> core::ops::IndexMut<(isize, isize, usize)> for CartesianDataFrame2D<S> {
    /// Exclusively borrows element at (i,j,n). The valid cells are indexed from zero 
    /// and ghost cells at the lower side of the domain are indexed with negative numbers
    /// This doesn't check that the index is valid as indices should only be used within
    /// iterators, which will only yield valid indices.
    fn index_mut(&mut self, indx: (isize,isize,usize)) -> &mut S {
        let (i,j,n) = indx;
        let p = get_flat_index_unchecked(i,j,n,&self.n_grown,self.n_ghost);
        &mut self.data[p]
    }
}


impl <S> core::ops::Index<(isize, isize, usize)> for CartesianDataFrame2D<S>{
    type Output = S;

    /// Borrows the element at (i,j,n). The valid cells are indexed from zero
    /// and ghost cells at the lower side of the domain are indexed with negative
    /// numbers.
    fn index(&self, indx: (isize, isize, usize) ) -> & Self::Output {
        let (i,j,n) = indx;
        let p = get_flat_index_unchecked(i, j, n, &self.n_grown, self.n_ghost);

        &self.data[p]
    }
}



trait GhostCells{
    /// checks if the cell at (i,j,k) contains a valid or a ghost cell.
    /// Returns true if valid, and returns false if ghost
    fn ij_is_valid_cell(&self, i: isize, j: isize) -> bool;

    /// check if the cell at p contains a valid or ghost cell. Returns same as ijk_is_valid_cell
    fn p_is_valid_cell(&self, p: usize) -> bool;
}


impl<S> GhostCells for CartesianDataFrame2D<S>
where
    S: Clone + Default
{
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


impl<S> CartesianDataFrame2D<S>
where
    S: Clone + Default + Copy + std::fmt::Debug
{
    fn collect_typed_cells(&mut self, side: Vec<CellType>) -> Vec<S> {
        self.iter()
            .ghost(side, IterComp::All)
            .cloned()
            .collect()
    }

    fn fill_prescribed_bc(&mut self, values: &Vec<S>, side: Vec<CellType>, comp: IterComp) {
        for (cell, value) in self.iter_mut().ghost(side, comp).zip(values.iter()){
            *cell = *value;
        }
    }
}

// Some functionality for Real CartesianDataFrames
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


impl BoundaryCondition<Real> for CartesianDataFrame2D<Real>{
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
                        unimplemented!("Reflect boundary condition not supported yet")
                    }
                    BcType::Prescribed(values) => {
                        // figure out the side we are working on
                        let mut sides: Vec<CellType> = Vec::new();
                        if i_dim == 0{
                            sides = vec![
                                CellType::Ghost(Loc::West),
                                CellType::Ghost(Loc::NorthWest),
                                CellType::Ghost(Loc::SouthWest)
                            ];
                        }
                        else if i_dim == 1 {
                            sides = vec![
                                CellType::Ghost(Loc::South),
                                CellType::Ghost(Loc::SouthWest),
                                CellType::Ghost(Loc::SouthEast)
                            ];
                        }

                        // fill the boundary condition
                        self.fill_prescribed_bc(values, sides, IterComp::Comps(vec![i_comp]));
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
                        unimplemented!("Reflect boundary condition not supported yet")
                    }
                    BcType::Prescribed(values) => {
                        let mut sides: Vec<CellType> = Vec::new();
                        if i_dim == 0 {
                            sides = vec![
                                CellType::Ghost(Loc::East),
                                CellType::Ghost(Loc::NorthEast),
                                CellType::Ghost(Loc::SouthEast)];
                        }
                        else if i_dim == 1 {
                            sides = vec![
                                CellType::Ghost(Loc::North),
                                CellType::Ghost(Loc::NorthEast),
                                CellType::Ghost(Loc::NorthWest)
                            ];
                        }
                        self.fill_prescribed_bc(values, sides, IterComp::Comps(vec![i_comp]));
                    }
                }
            }
        }
    }

}


/// Immutable iterator for `CartesianDataFrame`.
pub struct CartesianDataFrame2DIter<'a, S> {
    df: &'a CartesianDataFrame2D<S>,
    current_indx: usize,
}

fn cell_is_valid(cell_type: &CellType) -> bool {
    match cell_type {
        CellType::Valid(_) => { true }
        _ => { false }
    }
}

impl <'a, S> Iterator for CartesianDataFrame2DIter<'a, S>{
    type Item = &'a S;

    fn next(&mut self) -> Option<Self::Item> {
        // progress the current index to skip ghost cells
        while self.current_indx <= self.df.n_nodes - self.df.n_comp &&
            !cell_is_valid(&self.df.cell_type[self.current_indx])
            // !(match self.df.cell_type[self.current_indx]{CellType::Valid(_) => {true} _ => {false}})
        {
            self.current_indx += 1;
        }

        // check if the next item exists
        if self.current_indx <= self.df.n_nodes - self.df.n_comp &&
            cell_is_valid(&self.df.cell_type[self.current_indx])
            //(match self.df.cell_type[self.current_indx]{CellType::Valid(_) => {true} _ => {false}}){
                    //&& self.df.cell_type[self.current_indx] == CellType::Valid(_) {

        {
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
pub struct CartesianDataFrame2DIterMut<'a, S> {
    df: &'a mut CartesianDataFrame2D<S>,
    current_indx: usize,
}
impl <'a, S> Iterator for CartesianDataFrame2DIterMut<'a, S> {
    type Item = &'a mut S;

    // this is a safe function wrapping unsafe code. Rust cannot guarantee that it is safe, but
    // in practice it can be, and I'm pretty sure that it should be safe for everything we will 
    // ever want to do
    fn next(&mut self) -> Option<Self::Item> {
        // progress the current index to skip ghost cells
        while self.current_indx <= self.df.n_nodes - self.df.n_comp &&
            !(match self.df.cell_type[self.current_indx]{CellType::Valid(_) => {true} _ => {false}})
        {
            self.current_indx += 1;
        }

        // check if the next item exists
        if self.current_indx <= self.df.n_nodes-self.df.n_comp &&
            (match self.df.cell_type[self.current_indx]{CellType::Valid(_) => {true} _ => {false}})
        {
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


impl <'a, S> IntoIterator for &'a mut CartesianDataFrame2D<S>{
    type Item = &'a mut S;
    type IntoIter = CartesianDataFrame2DIterMut<'a, S>;

    fn into_iter(self) -> Self::IntoIter{
        Self::IntoIter {
            current_indx: 0,
            df: self,
        }
    }
}

impl<'a, S> IntoIterator for &'a CartesianDataFrame2D<S>{
    type Item = &'a S;
    type IntoIter = CartesianDataFrame2DIter<'a, S>;

    fn into_iter(self) -> Self::IntoIter{
        Self::IntoIter{
            current_indx: 0,
            df: self,
        }
    }
}


/// Struct to mutably iterate over a data frame, with the current index and the data itself
pub struct EnumerateIndex<'a, S>{
    iter: CartesianDataFrame2DIterMut<'a, S>,
}

/// Turns `CartesianDataFrameIter` into a mutable iterator which enumerates over the current index
pub trait IndexEnumerable <'a, S> {
    fn enumerate_index(self) -> EnumerateIndex<'a, S>;
}
impl <'a,S> IndexEnumerable <'a, S> for CartesianDataFrame2DIterMut<'a,S> {
    fn enumerate_index(self) -> EnumerateIndex<'a,S> {
        EnumerateIndex{
            iter: self,
        }
    }
}

impl <'a, S> Iterator for EnumerateIndex<'a,S>
    where S: Clone + Default
{
    type Item = ((isize, isize, usize), &'a mut S);

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
pub struct EnumeratePos<'a,S>{
    iter: CartesianDataFrame2DIterMut<'a,S>,
}

/// Turns `CartesianDataFrameIter` into a mutable iterator which enumerates the position 
/// of the data
pub trait PosEnumerable <'a,S> {
    fn enumerate_pos(self) -> EnumeratePos<'a,S>;
}

impl <'a,S> PosEnumerable <'a,S> for CartesianDataFrame2DIterMut<'a,S> {
    fn enumerate_pos(self) -> EnumeratePos<'a,S> {
        EnumeratePos{
            iter: self,
        }
    }
}

impl <'a, S> Iterator for EnumeratePos<'a,S>
    where S: Clone + Default
{
    type Item = ((Real, Real, usize), &'a mut S);

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

/// The components to iterate over
#[derive(Debug, Clone)]
pub enum IterComp{
    Comps(Vec<usize>),
    All,
}

impl IterComp {
    fn contains(&self, comp: usize) -> bool {
        match self {
            IterComp::Comps(cmps) => {cmps.contains(&comp)}
            IterComp::All => { true }
        }
    }
}

pub struct CartesianDataFrameGhostIter<'a, S> {
    df: &'a CartesianDataFrame2D<S>,
    current_indx: usize,
    side: Vec<CellType>,
    comps: IterComp,
}

/// mutable iterator over the east cells in `CartesianDataFrame`.
pub struct CartesianDataFrameGhostIterMut<'a,S> {
    df: &'a mut CartesianDataFrame2D<S>,
    current_indx: usize,
    side: Vec<CellType>,
    comps: IterComp,
}

impl <'a, S> Iterator for CartesianDataFrameGhostIterMut<'a,S>
    where S: Clone + Default
{
    type Item = &'a mut S;

    // this is a safe function wrapping unsafe code. Rust cannot guarantee that it is safe, but
    // in practice it can be, and I'm pretty sure that it should be safe for everything we will 
    // ever want to do
    fn next(&mut self) -> Option<Self::Item> {
        // progress the current index to next on the given side
        while self.current_indx < self.df.n_nodes &&
            // check if we want the side of the current cell
            (!(self.side.contains(&self.df.cell_type[self.current_indx])) ||
            // check if the component number is correct
            !(self.comps.contains(self.df.get_ij(self.current_indx).unwrap().2)))
        {
            self.current_indx += 1;
        }

        // check if the next item exists
        if self.current_indx < self.df.n_nodes {
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

impl <'a, S> Iterator for CartesianDataFrameGhostIter<'a, S>
    where S: Clone + Default
{
    type Item = &'a S;

    fn next(&mut self) -> Option<Self::Item> {
        // progress the current index to next on the given side
        while self.current_indx < self.df.n_nodes &&
            (!(self.side.contains(&self.df.cell_type[self.current_indx])) ||
            !(self.comps.contains(self.df.get_ij(self.current_indx).unwrap().2)))
        {
            self.current_indx += 1;
        }

        // check if the next item exists
        if self.current_indx < self.df.n_nodes{
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

pub trait GhostIteratorMut <'a,S> {
    fn ghost(self, side: Vec<CellType>, comps: IterComp) -> CartesianDataFrameGhostIterMut<'a,S>;
}

pub trait GhostIterator <'a, S> {
    fn ghost(self, side: Vec<CellType>, comps: IterComp) -> CartesianDataFrameGhostIter<'a, S>;
}

impl <'a,S> GhostIteratorMut <'a,S> for CartesianDataFrame2DIterMut<'a,S> {
    fn ghost(self, side: Vec<CellType>, comps: IterComp) -> CartesianDataFrameGhostIterMut<'a,S> {
        CartesianDataFrameGhostIterMut{
            current_indx: 0,
            df: self.df,
            side,
            comps,
        }   
    }
}

impl <'a, S> GhostIterator <'a, S> for CartesianDataFrame2DIter<'a, S> {
    fn ghost(self, side: Vec<CellType>, comps: IterComp) -> CartesianDataFrameGhostIter<'a, S> {
        CartesianDataFrameGhostIter{
            current_indx: 0,
            df: self.df,
            side,
            comps,
        }
    }
}


/// Struct to iterate over a data frame, with the current index and the data itself
pub struct EnumerateGhostIndex<'a,S>{
    iter: CartesianDataFrameGhostIterMut<'a,S>,
}

/// Turns `CartesianDataFrameIter` into an iterator which enumerates over the current index
pub trait GhostIndexEnumerable <'a,S> {
    fn enumerate_index(self) -> EnumerateGhostIndex<'a,S>;
}

impl <'a,S> GhostIndexEnumerable <'a,S> for CartesianDataFrameGhostIterMut<'a,S> {
    fn enumerate_index(self) -> EnumerateGhostIndex<'a,S> {
        EnumerateGhostIndex{
            iter: self,
        }
    }
}

impl <'a, S> Iterator for EnumerateGhostIndex<'a,S>
    where S: Clone + Default
{
    type Item = ((isize, isize, usize), &'a mut S);

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
    fn test_iter_comp(){
        let cmps = IterComp::Comps(vec![1,2]);
        assert!(cmps.contains(1) == true);
        assert!(cmps.contains(5) == false);
    }

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
            assert_eq!(df2_iter.next(), Some(& 2.0));
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
                CellType::Ghost(Loc::SouthWest),
                CellType::Ghost(Loc::South),
                CellType::Ghost(Loc::South),
                CellType::Ghost(Loc::South),
                CellType::Ghost(Loc::SouthEast),
                CellType::Ghost(Loc::West),
                CellType::Valid(Loc::SouthWest),
                CellType::Valid(Loc::South),
                CellType::Valid(Loc::SouthEast),
                CellType::Ghost(Loc::East),
                CellType::Ghost(Loc::West),
                CellType::Valid(Loc::West),
                CellType::Valid(Loc::Regular),
                CellType::Valid(Loc::East),
                CellType::Ghost(Loc::East),
                CellType::Ghost(Loc::West),
                CellType::Valid(Loc::NorthWest),
                CellType::Valid(Loc::North),
                CellType::Valid(Loc::NorthEast),
                CellType::Ghost(Loc::East),
                CellType::Ghost(Loc::NorthWest),
                CellType::Ghost(Loc::North),
                CellType::Ghost(Loc::North),
                CellType::Ghost(Loc::North),
                CellType::Ghost(Loc::NorthEast)
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
    fn test_prescribed_two_ghost() {
        let u1 = CartesianMesh2D::new(RealVec2([0.0, 0.0]), RealVec2([6.0, 6.0]), UIntVec2([3,3]));
        let values = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0];
        let bc = BCs::new(
            vec![
                ComponentBCs::new(
                    vec![BcType::Prescribed(values.clone()), BcType::Prescribed(values.clone())],
                    vec![BcType::Prescribed(values.clone()), BcType::Prescribed(values.clone())],
                )
            ]
        );


        let mut cdf = CartesianDataFrame2D::new_from(&u1, bc, 1, 2);
        cdf.fill_bc();
        assert_eq!(
            cdf.data,
            vec![1.0,  2.0,  3.0,  4.0,  5.0,  6.0,  7.0,
                 8.0,  9.0, 10.0, 11.0, 12.0, 13.0, 14.0,
                 5.0,  6.0,  0.0,  0.0,  0.0,  5.0,  6.0,
                 7.0,  8.0,  0.0,  0.0,  0.0,  7.0,  8.0,
                 9.0, 10.0,  0.0,  0.0,  0.0,  9.0, 10.0,
                 1.0,  2.0,  3.0,  4.0,  5.0,  6.0,  7.0,
                 8.0,  9.0, 10.0, 11.0, 12.0, 13.0, 14.0]
        )
    }

    #[test]
    fn test_ghost_iter_mut_two_ghost () {
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
        //     0.0,  0.0, 6.0, 4.0, 2.0, 0.0,  0.0,
        //     0.0,  0.0, 6.0, 4.0, 2.0, 0.0,  0.0,
        //    -4.0, -2.0, 2.0, 4.0, 6.0, 8.0, 10.0,
        //    -4.0, -2.0, 2.0, 4.0, 6.0, 8.0, 10.0,
        //    -4.0, -2.0, 2.0, 4.0, 6.0, 8.0, 10.0,
        //     0.0,  0.0, 4.0, 2.0, 0.0, 0.0,  0.0,
        //     0.0,  0.0, 4.0, 2.0, 0.0, 0.0,  0.0
        //]
        let mut east_iter = cdf.iter_mut().ghost(vec![CellType::Ghost(Loc::East)], IterComp::All);
        assert_eq!(east_iter.next(), Some(&mut 8.0));
        assert_eq!(east_iter.next(), Some(&mut 10.0));
        assert_eq!(east_iter.next(), Some(&mut 8.0));
        assert_eq!(east_iter.next(), Some(&mut 10.0));
        assert_eq!(east_iter.next(), Some(&mut 8.0));
        assert_eq!(east_iter.next(), Some(&mut 10.0));
        assert_eq!(east_iter.next(), None);
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
        //     0.0,  0.0, 6.0, 4.0, 2.0, 0.0,  0.0,
        //     0.0,  0.0, 6.0, 4.0, 2.0, 0.0,  0.0,
        //    -4.0, -2.0, 2.0, 4.0, 6.0, 8.0, 10.0,
        //    -4.0, -2.0, 2.0, 4.0, 6.0, 8.0, 10.0,
        //    -4.0, -2.0, 2.0, 4.0, 6.0, 8.0, 10.0,
        //     0.0,  0.0, 4.0, 2.0, 0.0, 0.0,  0.0,
        //     0.0,  0.0, 4.0, 2.0, 0.0, 0.0,  0.0
        //]
        //


        let mut east_iter = cdf.iter().ghost(vec![CellType::Ghost(Loc::East)], IterComp::All);
        assert_eq!(east_iter.next(), Some(& 8.0));
        assert_eq!(east_iter.next(), Some(& 10.0));
        assert_eq!(east_iter.next(), Some(& 8.0));
        assert_eq!(east_iter.next(), Some(& 10.0));
        assert_eq!(east_iter.next(), Some(& 8.0));
        assert_eq!(east_iter.next(), Some(& 10.0));
        assert_eq!(east_iter.next(), None);

        let mut south_iter = cdf.iter().ghost(vec![CellType::Ghost(Loc::South)], IterComp::All);
        assert_eq!(south_iter.next(), Some(& 4.0));
        assert_eq!(south_iter.next(), Some(& 2.0));
        assert_eq!(south_iter.next(), Some(& 0.0));
        assert_eq!(south_iter.next(), Some(& 4.0));
        assert_eq!(south_iter.next(), Some(& 2.0));
        assert_eq!(south_iter.next(), Some(& 0.0));
        assert_eq!(south_iter.next(), None);
    }

    #[test]
    fn test_ghost_collect_two_edge() {
        let u1 = CartesianMesh2D::new(
            RealVec2([0.0, 0.0]),
            RealVec2([10.0, 10.0]),
            UIntVec2([5, 5])
        );

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

        let north = cdf.collect_typed_cells(
            vec![
                CellType::Valid(Loc::North),
                CellType::Valid(Loc::NorthEast),
                CellType::Valid(Loc::NorthWest)
            ]
        );
        assert_eq!(north, vec![2.0, 4.0, 6.0, 8.0, 10.0, 2.0, 4.0, 6.0, 8.0, 10.0]);

        let south = cdf.collect_typed_cells(
            vec![
                CellType::Valid(Loc::South),
                CellType::Valid(Loc::SouthWest),
                CellType::Valid(Loc::SouthEast)
            ]
        );
        assert_eq!(south, vec![2.0, 4.0, 6.0, 8.0, 10.0, 2.0, 4.0, 6.0, 8.0, 10.0]);

        let east = cdf.collect_typed_cells(
            vec![
                CellType::Valid(Loc::East),
                CellType::Valid(Loc::NorthEast),
                CellType::Valid(Loc::SouthEast)
            ]
        );
        assert_eq!(east, vec![8.0, 10.0, 8.0, 10.0, 8.0, 10.0, 8.0, 10.0, 8.0, 10.0]);

        let west = cdf.collect_typed_cells(
            vec![
                CellType::Valid(Loc::West),
                CellType::Valid(Loc::NorthWest),
                CellType::Valid(Loc::SouthWest)
            ]
        );
        assert_eq!(west, vec![2.0, 4.0, 2.0, 4.0, 2.0, 4.0, 2.0, 4.0, 2.0, 4.0]);
    }
}

