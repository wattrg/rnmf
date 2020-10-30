
/// Structure containing data to define a CartesianMesh
#[derive(Debug)]
pub struct CartesianMesh
{
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



impl CartesianMesh
{


    /// Generate new CartesianMesh from the lo and high corner, and the number of nodes
    pub fn new(lo: Vec<f64>, hi: Vec<f64>, n:Vec<usize>, dim: u32) -> CartesianMesh
    {
        let mut u = CartesianMesh
        {
            lo,
            hi,
            n,
            dx: vec![0.0; dim as usize],
            node_pos: Vec::new(),
            dim
        };

        for i_dim in 0..dim
        {
            for i in 0..u.n[i_dim as usize]
            {
                    u.dx[i_dim as usize] = (u.hi[i_dim as usize] - u.lo[i_dim as usize])/(u.n[i_dim as usize] as f64);
                    u.node_pos.push(u.lo[i_dim as usize] + 0.5*u.dx[i_dim as usize] + (i as f64)*u.dx[i_dim as usize]);
            }
        }
        u
    }
}

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
