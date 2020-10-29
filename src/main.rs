

mod mesh;
mod data_frame;

use data_frame::BCType;

/// function describing the initial condition
#[allow(unused_variables)]
fn initial_condition (x: f64) -> f64
{
    2.0
}

fn main() {
    println!("Hello, world!");

    let u1 = mesh::CartesianMesh::new(0.0, 10.0, 10);
    let mut cdf = data_frame::CartesianDataFrame::new_from(u1, 2); 
    //note df now owns u1, and u1 is invalid

    cdf.fill_ic(&initial_condition);
    cdf.fill_bc(&BCType::Dirichlet(1.0), &BCType::Dirichlet(1.0));
    println!("CartesianDataFrame = {:?}", cdf);
}
// Hello World!! Kind Regards Gabby <3 
