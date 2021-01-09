mod poisson;
mod model;

use rnmf::*;
use boundary_conditions::*;
use mesh::cartesian::*;
use std::env;
use poisson::*;
use io::{RnmfProgressBar};
//use rand::Rng;

fn main() {
    println!("Hello from the magnetistatics example!");

    // read command line argument
    let args: Vec<String> = env::args().collect();
    match args.len() > 1 {
        true => {}
        false => panic!("please provide the lua configuration file")
    }

    // read lua file
    let conf = config::read_lua(&args[1]).unwrap();
    let model_conf = model::read_lua(&args[1]).unwrap();



    // create the mesh 
    let mesh = CartesianMesh::new(
        vec![0.0, 0.0],                                     // lo corner
        vec![model_conf.length[0], model_conf.length[1]],   // hi corner
        vec![model_conf.n_cells[0], model_conf.n_cells[1]], // number of cells
        conf.dim                                            // dimensions
    );


    let bc = BCs::new(vec![
        ComponentBCs::new(
            // lo bc
            vec![BCType::Neumann(-model_conf.h_far[0]/mesh.dx[0]), 
                 BCType::Neumann(-model_conf.h_far[1]/mesh.dx[1])], 
            // hi bc
            vec![BCType::Neumann(-model_conf.h_far[0]/mesh.dx[0]), 
                 BCType::Neumann(-model_conf.h_far[1]/mesh.dx[1])]  
        )
    ]);

    // create the data frame
    let mut psi = CartesianDataFrame::new_from(&mesh, 1, 1);
    psi.fill_ic(|x,y,_,_| -> f64 {
        -model_conf.h_far[0]/mesh.dx[0]*x - model_conf.h_far[1]/mesh.dx[1]*y
        
    });
    psi.fill_bc(&bc);
    io::write_csv(&"./examples/magnetostatics/psi_initial", psi.clone());

    // begin solving
    println!("Overall progress:");
    let mut progress_bar = io::ProgressBar::create(model_conf.n_iter);
    let mut source = psi.clone();
    for _ in 0..model_conf.n_iter{
        source = get_source(&psi, &model_conf);
        let psi_star = solve_poisson(&mut psi, &source, &model_conf, &bc);
        let diff = &psi_star + &(-1.0 * &psi);
        psi = &psi + &(model_conf.relax * diff);
        psi.fill_bc(&bc);
        progress_bar.increment(1);
    }
    progress_bar.finish();

    io::write_csv(&"./examples/magnetostatics/psi_final", psi);
    io::write_csv(&"./examples/magnetostatics/source", source);
    println!("Done.");

}