use rnmf::*;
use boundary_conditions::*;
use rnmf::mesh::cartesian2d::{
    mesh::CartesianMesh2D, 
    dataframe::{
        CartesianDataFrame2D, 
        BoundaryCondition2D
    }
};
use std::env;

fn main()  {
    // read command line argument
    let args: Vec<String> = env::args().collect();
    match args.len() > 1 {
        true => {}
        false => panic!("please provide the lua configuration file")
    }

    // define the boundary conditions
    let bc = BCs::new(vec![
        ComponentBCs::new(
                // x direction
                vec![BcType::Dirichlet(0.0), BcType::Dirichlet(1.0)], // lo

                // y direction
                vec![BcType::Dirichlet(0.0), BcType::Dirichlet(1.0)], // hi
        )]
    );

    // create the mesh 
    let u1 = CartesianMesh2D::new(RealVec2([0.0, 0.0]), RealVec2([6.0, 6.0]), UIntVec2([3, 3]));
   
    // create some variables on top of the mesh
    let mut variable1 = CartesianDataFrame2D::new_from(&u1, bc.clone(), 2, 1);
    let mut variable2 = CartesianDataFrame2D::new_from(&u1, bc, 1, 1);

    // fill the initial conditions
    variable1.fill_ic(|x,y,n| x * y * (n + 1) as Real);
    variable2.fill_ic(|x,y,_| x + y);


    variable2.fill_bc();

    // print out the variables
    println!("variable 1 = {:?}", variable1);
    println!();
    println!("variable 2 = {:?}", variable2);



}

// Hello World!! Kind Regards Gabby <3

