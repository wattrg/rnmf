
/// Structure containing data to define a CartesianMesh
#[derive(Debug)]
pub struct CartesianMesh
{
    /// The position of the nodes in the mesh
    pub node_pos: Vec<f64>,

    /// The lower corner of the mesh
    pub lo: f64,

    /// The upper corner of the mesh
    pub hi: f64,

    /// The number of cells in nodes in each direction
    pub n : usize,

    /// The distance between each node in each
    pub dx: f64
}



impl CartesianMesh
{


    /// Generate new CartesianMesh from the lo and high corner, and the number of nodes
    pub fn new(lo: f64, hi: f64, n:usize) -> CartesianMesh
    {
        let mut u = CartesianMesh
        {
            lo,
            hi,
            n,
            dx: (hi-lo)/(n as f64),
            node_pos: Vec::new(),
        };

        for i in 0..n
        {
            u.node_pos.push(u.lo + 0.5*u.dx + (i as f64)*u.dx)
        }
        u
    }
}
