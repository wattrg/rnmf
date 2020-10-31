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
    pub dim: u32
}

impl CartesianMesh {
    /// Generate new CartesianMesh from the lo and high corner, and the number of nodes
    pub fn new(lo: Vec<f64>, hi: Vec<f64>, n:Vec<usize>, dim: u32) -> CartesianMesh
    {
        // This will have to be updated to calculate the node positions for meshes in 1, 2, and 3D eventually,
        // but its time to study XD
        let mut dx = vec![0.0; dim as usize];
        for i_dim in 0..dim
        {
            //for i in 0..n[i_dim as usize]
            //{
                    dx[i_dim as usize] = (hi[i_dim as usize] - lo[i_dim as usize])/(n[i_dim as usize] as f64);
                    //u.node_pos.push(u.lo[i_dim as usize] + 0.5*u.dx[i_dim as usize] + (i as f64)*u.dx[i_dim as usize]);
            //}
        }
        let node_pos: Vec<f64> = (0..n[0]).map(|i| lo[0] + ((i as f64) + 0.5) * dx[0]).collect();
        CartesianMesh
        {
            lo,
            hi,
            n,
            dx,
            node_pos,
            dim
        }
    }
}

// unused for the moment, but will form the basis for indexing the mesh
// Will need to make a iterators to iterate over the data
#[allow(dead_code)]
impl CartesianMesh
{
    /// Retrieves the element at (i,j,k,n)
    fn index_mut(&mut self, i: usize, j: usize, k: usize) -> f64
    {
        self.node_pos[
            i + 
            self.n[1]*j +
            self.n[2]*self.n[1]*k
        ]
    }
}
