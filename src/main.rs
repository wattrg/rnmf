use rnmf::*;
use boundary_conditions::*;
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
    let mut variable2 = CartesianDataFrame::new_from(&u1, 1, 1);

    // fill the initial conditions
    variable1.fill_ic(|x,y,_,n| x * y * (n + 1) as Real);
    variable2.fill_ic(|x,y,_,_| x + y);

    let bc = BCs::new(vec![
        ComponentBCs::new(
                // x direction
                vec![BCType::Dirichlet(0.0), BCType::Dirichlet(1.0)], // lo

                // y direction
                vec![BCType::Dirichlet(0.0), BCType::Dirichlet(1.0)], // hi
        )]
    );

    variable2.fill_bc(&bc);

    // print out the variables
    println!("variable 1 = {:?}", variable1);
    println!();
    println!("variable 2 = {:?}", variable2);



}

// Hello World!! Kind Regards Gabby <3

