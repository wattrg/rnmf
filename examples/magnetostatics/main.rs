mod poisson;
mod model;

use rnmf::*;
use boundary_conditions::*;
use mesh::cartesian2d::*;
use std::env;
use poisson::*;
use io::*;
use colored::*;
use crate::config::{UserConfig};

fn main() {
    // read command line argument
    let args: Vec<String> = env::args().collect();

    // generate instance of the user model
    let model_instance = model::UserModel::new();

    // configure the simulation
    let (conf, out_cb) = config::init(args, model_instance, model::output_callbacks).unwrap();

    // create the mesh 
    let mesh = CartesianMesh2D::new(
        [0.0, 0.0],                                   // lo corner
        [conf.geom.length[0], conf.geom.length[1]],   // hi corner
        [conf.geom.n_cells[0], conf.geom.n_cells[1]]  // number of cells
    );


    // create the boundary conditions
    let bc = BCs::new(vec![
        ComponentBCs::new(
            // lo bc
            vec![BCType::Neumann(-conf.model.h_far[0]/mesh.dx[0]), 
                 BCType::Neumann(-conf.model.h_far[1]/mesh.dx[1])], 
            // hi bc
            vec![BCType::Neumann(-conf.model.h_far[0]/mesh.dx[0]), 
                 BCType::Neumann(-conf.model.h_far[1]/mesh.dx[1])]  
        )
    ]);

    // create the data frame
    let mut psi = CartesianDataFrame2D::new_from(&mesh, 1, 1);
    psi.fill_ic(|x,y,_| -> Real {
        -conf.model.h_far[0]/mesh.dx[0]*x - conf.model.h_far[1]/mesh.dx[1]*y
        
    });
    psi.fill_bc(&bc);
    
    println!("Calculating solution:");
    let mut progress_bar = io::ProgressBar::create(conf.model.n_iter);

    io::write_vtk(&"./examples/magnetostatics/", &"ferro_droplet", 0, &psi, &out_cb, &progress_bar);

    // begin solving
    let mut sum_diff: Real = 1.0;
    let mut iter = 0;
    while iter < conf.model.n_iter && (sum_diff > conf.model.tol){
        let source = get_source(&psi, &conf);
        let psi_star = solve_poisson(&mut psi, &source, &conf, &bc);
        let diff = &psi_star + &(-1.0 * &psi);
        sum_diff = sum(&diff);
        psi = &psi + &(conf.model.relax * diff);
        psi.fill_bc(&bc);
        progress_bar.increment(1);
        iter += 1;
    }
    progress_bar.finish();

    io::write_vtk(&"./examples/magnetostatics/", &"ferro_droplet", conf.model.n_iter, &psi, &out_cb, &progress_bar);
    println!("sum diff = {}", sum_diff);
    println!("{}", "Done.".green().bold());

}