mod data_frame;
mod mesh;

use data_frame::{BCType, BoundaryCondition};

fn main() {
    println!("Hello, world!");

    let u1 = mesh::CartesianMesh::new(0.0, 10.0, 10);
    let mut variable1 = data_frame::CartesianDataFrame::new_from(&u1, 2);
    let mut variable2 = data_frame::CartesianDataFrame::new_from(&u1, 2);

    variable1.fill_ic(|_| 2.0);
    variable1.fill_bc(BCType::Dirichlet(1.0), BCType::Dirichlet(1.0));

    variable2.fill_ic(|_| 2.0);
    variable2.fill_bc(BCType::Dirichlet(0.0), BCType::Dirichlet(0.0));

    println!("Variable 1 = {:?}", variable1);
    println!("Variable 2 = {:?}", variable2);
}
// Hello World!! Kind Regards Gabby <3
