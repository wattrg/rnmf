/// ## Module containing mesh implementations
/// Each type of mesh is a sub-module within this module. The currently available mesh types are
/// * Cartesian
mod mesh;

/// ## Configuration module
/// Contains simulation wide settings which define the simulation
mod config;

/// ## Module defining the available boundary conditions
/// It defines the `BoundaryCondition` trait, which all DataFrames must have to fill the boundary
/// conditions. It also contains an enum `BCType` which lists the available boundary conditions.
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

    println!("{:?}", variable1);
    {
        let variable1_iter = CartesianDataFrameIter{
            current_indx: 0,
            df: &mut variable1,
        };

        for d in variable1_iter {
            *d += 1.0;
        }
    }

    println!("{:?}", variable1)


}
// Hello World!! Kind Regards Gabby <3

