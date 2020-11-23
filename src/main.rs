mod mesh;
mod config;
mod boundary_conditions;

use boundary_conditions::{BCType, BoundaryCondition};
use crate::mesh::cartesian::*;


fn main() {
    println!("Hello, world!");

    let config = config::Config{
        dim: 1 as usize,
    };

    let u1 = CartesianMesh::new(vec![0.0, 0.0], vec![10.0, 10.0], vec![10, 10], config.dim);
    let mut variable1 = CartesianDataFrame::new_from(&u1, 1, 2);
    let mut variable2 = CartesianDataFrame::new_from(&u1, 1, 2);

    variable1.fill_ic(|_,_,_| 2.0);
    variable1.fill_bc(BCType::Dirichlet(1.0), BCType::Dirichlet(1.0));

    variable2.fill_ic(|_,_,_| 2.0);
    variable2.fill_bc(BCType::Dirichlet(0.0), BCType::Dirichlet(0.0));

    println!("node positions = {:?}", u1);

    println!("Variable 1 = {:?}", variable1);
    println!("Variable 2 = {:?}", variable2);

}
// Hello World!! Kind Regards Gabby <3

