

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


impl CartesianMesh {
    /// Generate new CartesianMesh from the lo and high corner, 
    /// and the number of nodes
    pub fn new(lo: Vec<f64>, hi: Vec<f64>, n:Vec<usize>, dim: usize) 
        -> CartesianMesh
    {
        let mut dx = vec![0.0; dim as usize];
        for i_dim in 0..dim as usize{
            dx[i_dim] = (hi[i_dim] - lo[i_dim])/(n[i_dim] as f64);
        }

        let mut n_nodes = n[0];
        if dim >= 2{
            n_nodes *= n[1];

            if dim == 3 {
                n_nodes *= n[2];
            }
        }

        n_nodes *= dim;
        
        // let node_pos: Vec<f64> = (0..n[0])
        //                             .map(|i| lo[0] + ((i as f64) + 0.5) * dx[0])
        //                             .collect();


        let node_pos = Vec::with_capacity(n_nodes);


        

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
        cm.node_pos = (0..n_nodes).map(
            |p: usize| -> f64 {
                let (i,j,k,n) = cm.get_ijk(p).unwrap();
                match dim {
                    1 => {
                        cm.lo[n] + ((i as f64) + 0.5) * cm.dx[n]
                    }
                    2 => {
                        cm.lo[n] + ((j as f64) + 0.5) * cm.dx[n]
                    }
                    3 => {
                        cm.lo[n] + ((k as f64) + 0.5) * cm.dx[n]
                    }
                    _ => {
                        panic!("{}D not supported!");
                    }
                }
                
            }).collect();
        cm

    }
}

/// Macro which you can pass the mesh and an index and it will return the flattened index
#[macro_export]
macro_rules! index {
    ($self:ident, $i:expr ) => {
        $i
    };
    
    ($self:ident, $i:expr, $j:expr) => {
        $i + $self.n[0] * $j
    };
    
    ($self:ident, $i:expr, $j:expr, $k:expr) => {
        $i * $self.n[0] * $j * $self.n[1] * $self.n[2] * $k
    };
}

// unused for the moment, but will form the basis for indexing the mesh
impl CartesianMesh
{
    /// Retrieves the element at (i,j,k,n)
    #[allow(dead_code)]
    fn index_mut(&mut self, i: usize, j: usize, k: usize) -> f64
    {
        let mut p = i;
        if self.dim >= 2 {
            p += self.n[1]*j;
        
            if self.dim == 3 {
                p += self.n[2]*self.n[1]*k;
            }
        }
        self.node_pos[p]
    }

    /// Returns the un-flattened index
    #[allow(dead_code)]
    fn get_ijk(& self, p: usize) -> Option<(usize, usize, usize, usize)>{
        if p >= self.n_nodes {
            None
        }
        else{
            // the idea with the closures is that I'll implement lazy evaluation for each 
            // closure following https://doc.rust-lang.org/book/ch13-01-closures.html, so that we
            // aren't performing any unnecessary or repeated calculations
            let d = || {
                p % self.dim
            };
            
            let i = || {
                ((p-d())/self.dim) % self.n[0]
            };
            let j = || {
                ((p - d() - i()*self.dim)/(self.dim * self.n[0])) % self.n[1]
            };
            let k = || {
                (p- j()*self.n[0]*self.dim - i()*self.dim - d())/(self.dim * self.n[0] * self.n[1])
            };
            

            match self.dim{
                1 => {
                    Some((i(), 0, 0, d()))
                }
                2 => {
                    Some((i(), j(), 0, d()))
                }
                3 => {
                    Some((i(), j(), k(), d()))
                }
                _ => {
                    panic!("{}D not supported!", self.dim);
                }
            }
        }
    }
}


/// Immutable Iterator struct for Cartesian
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
        if self.current_indx >= self.mesh.node_pos.len() {
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




#[cfg(test)]
mod tests{
    use super::*;

    #[test]
    fn test_macro (){
        let u1 = CartesianMesh::new(vec![0.0, 0.0], vec![10.0, 10.0], vec![5, 2], 2);
        assert_eq!(index!(u1, 2, 3), 17);
    }

    #[test]
    fn test_iterator () {
        let m = CartesianMesh::new(vec![0.0], vec![10.0], vec![5], 1);
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
        let m1 = CartesianMesh::new(vec![0.0], vec![10.0], vec![5], 1);
        let m2 = CartesianMesh::new(vec![0.0,0.0], vec![10.0,8.0], vec![5,4], 2);
        let m3 = CartesianMesh::new(vec![0.0,0.0,0.0], vec![10.0,8.0,10.0], vec![5,4,5], 3);

        assert_eq!(m1.get_ijk(3).unwrap(), (3,0,0,0));
        assert_eq!(m1.get_ijk(5), Option::None);

        assert_eq!(m2.get_ijk(3).unwrap(), (1,0,0,1));
        assert_eq!(m2.get_ijk(15).unwrap(), (2,1,0,1));
        assert_eq!(m2.get_ijk(40), Option::None);

        assert_eq!(m3.get_ijk(3).unwrap(), (1,0,0,0));
        assert_eq!(m3.get_ijk(22).unwrap(), (2,1,0,1));
        assert_eq!(m3.get_ijk(60).unwrap(), (0,0,1,0));
        assert_eq!(m3.get_ijk(300), Option::None);
    }
}

