

mod mesh;
mod data_frame;

use data_frame::{BCType, BoundaryCondition};

/// function describing the initial condition
#[allow(unused_variables)]
fn initial_condition (x: f64, y: f64, z: f64) -> f64
{
    2.0
}

fn main() {
    println!("Hello, world!");

    let u1 = mesh::CartesianMesh::new(vec![0.0], vec![10.0], vec![10], 1);
    let mut variable1 = data_frame::CartesianDataFrame::new_from(&u1, 1, 2);
    //let mut variable2 = data_frame::CartesianDataFrame::new_from(&u1, 1, 2);

    variable1.fill_ic(&initial_condition);
    variable1.fill_bc(&BCType::Dirichlet(1.0), &BCType::Dirichlet(1.0));

    //variable2.fill_ic(&initial_condition);
    //variable2.fill_bc(&BCType::Dirichlet(0.0), &BCType::Dirichlet(0.0));

    println!("Variable 1 = {:?}", variable1);
    //println!("Variable 2 = {:?}", variable2);
}
// Hello World!! Kind Regards Gabby <3 
