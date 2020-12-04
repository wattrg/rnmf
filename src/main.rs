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

#[allow(unused_imports)]
use boundary_conditions::{BCType, BoundaryCondition};
use crate::mesh::cartesian::*;
use std::env;


fn main()  {
    // read command line argument
    let args: Vec<String> = env::args().collect();
    match args.len() > 1 {
        true => {}
        false => panic!("please provide the lua configuration file")
    }
    let conf = config::read_lua(&args[1]).unwrap();

    // create the mesh 
    let u1 = CartesianMesh::new(vec![0.0, 0.0], vec![6.0, 6.0], vec![3, 3], conf.dim);
   
    // create some variables on top of the mesh
    let mut variable1 = CartesianDataFrame::new_from(&u1, 2, 1);
    let mut variable2 = CartesianDataFrame::new_from(&u1, 1, 2);

    // fill the initial conditions
    variable1.fill_ic(|x,y,_,n| x * y * (n + 1) as f64);
    variable2.fill_ic(|x,y,_,_| x + y);

    // print out the variables
    println!("variable 1 = {:?}", variable1);
    println!("");
    println!("variable 2 = {:?}", variable2);



}

// Hello World!! Kind Regards Gabby <3

