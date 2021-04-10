
use super::*;
use std::rc::Rc;


/// Structure containing data to define a CartesianMesh
#[derive(Debug, Clone, Default)]
pub struct CartesianMesh2D {
    /// The position of the nodes in the mesh
    pub node_pos: Vec<Real>,

    /// The lower corner of the mesh
    pub lo: RealVec2,

    /// The upper corner of the mesh
    pub hi: RealVec2,

    /// The number of cells in nodes in each direction
    pub n : UIntVec2,

    /// The distance between each node in each
    pub dx: RealVec2,

    /// The number of dimensions
    pub dim: usize,
    
    /// The number of nodes in the mesh
    n_nodes: usize
}



impl CartesianMesh2D {
    /// Generate new CartesianMesh from the lo and high corner, 
    /// and the number of nodes
    pub fn new(lo: RealVec2, hi: RealVec2, n: UIntVec2) -> Rc<CartesianMesh2D>
    {
        let mut dx = [0.0; 2];
        for i_dim in 0..2_usize{
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
            dx: RealVec2(dx),
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




#[cfg(test)]
mod tests{
    use super::*;

    #[test]
    fn mesh_iterator () {

        let m2 = CartesianMesh2D::new(
            RealVec2([0.0, 0.0]), RealVec2([6.0, 6.0]), UIntVec2([3, 3])
        );
        let mut m2_iter = m2.into_iter();
        assert_eq!(m2_iter.next().unwrap(), vec![& 1.0, & 1.0]);
        assert_eq!(m2_iter.next().unwrap(), vec![& 3.0, & 1.0]);
        assert_eq!(m2_iter.next().unwrap(), vec![& 5.0, & 1.0]);
        assert_eq!(m2_iter.next().unwrap(), vec![& 1.0, & 3.0]);
        assert_eq!(m2_iter.next().unwrap(), vec![& 3.0, & 3.0]);
        assert_eq!(m2_iter.next().unwrap(), vec![& 5.0, & 3.0]);
        assert_eq!(m2_iter.next().unwrap(), vec![& 1.0, & 5.0]);
        assert_eq!(m2_iter.next().unwrap(), vec![& 3.0, & 5.0]);
        assert_eq!(m2_iter.next().unwrap(), vec![& 5.0, & 5.0]);
        assert_eq!(m2_iter.next(), Option::None);
    }

    #[test]
    fn mesh_indexing () {

        let m2 = CartesianMesh2D::new(RealVec2([0.0,0.0]), RealVec2([6.0,6.0]), UIntVec2([3,3]));
        assert_eq!(m2.get_ij(16).unwrap(), (1,2,1));
        assert_eq!(m2.get_ij(40), Option::None);
        assert_eq!(m2.index(2,1,1), &3.0);

    }

    #[test]
    fn node_pos () {

        let m2 = CartesianMesh2D::new(
            RealVec2([0.0, 0.0]), RealVec2([10.0, 10.0]), UIntVec2([5,5])
        );
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
}