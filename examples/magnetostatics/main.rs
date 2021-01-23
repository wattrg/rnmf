mod poisson;
mod model;

use rnmf::*;
use boundary_conditions::*;
use mesh::cartesian2d::*;
use std::env;
use poisson::*;
use io::*;
use colored::*;
use std::process::Command;
use std::collections::HashMap;

fn main() {
    let version = env!("CARGO_PKG_VERSION");
    let commit = String::from_utf8(Command::new("git").args(&["rev-parse", "HEAD"]).output().unwrap().stdout).unwrap();

    print!("{}", format!("rnmf v {}, commit {}", version, commit).green().bold());
    println!("Hello from the magnetistatics example!");


    if !cfg!(feature = "disable_double") {
        println!("using double precision");
    }

    // read command line argument
    let args: Vec<String> = env::args().collect();
    match args.len() > 1 {
        true => {}
        false => panic!("please provide the lua configuration file")
    }

    // read lua file
    let model_conf = model::read_lua(&args[1]).unwrap();


    //let output: std::collections::HashMap<String, io::OutputCallBack> = 
    //[("psi", model::get_psi), ("h", model::get_h)].iter().cloned().collect();

    let mut output: HashMap<String, io::OutputCallBack> = HashMap::new();
    output.insert("psi".to_string(), model::get_psi);
    output.insert("h".to_string(), model::get_h);

    // create the mesh 
    let mesh = CartesianMesh2D::new(
        [0.0, 0.0],                                     // lo corner
        [model_conf.length[0], model_conf.length[1]],   // hi corner
        [model_conf.n_cells[0], model_conf.n_cells[1]]  // number of cells
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
    let mut psi = CartesianDataFrame2D::new_from(&mesh, 1, 1);
    psi.fill_ic(|x,y,_| -> Real {
        -model_conf.h_far[0]/mesh.dx[0]*x - model_conf.h_far[1]/mesh.dx[1]*y
        
    });
    psi.fill_bc(&bc);
    
    println!("Calculating solution:");
    let mut progress_bar = io::ProgressBar::create(model_conf.n_iter);

    io::write_vtk(&"./examples/magnetostatics/", &"ferro_droplet", 0, &psi, &output, &progress_bar);

    // begin solving
    for _ in 0..model_conf.n_iter{
        let source = get_source(&psi, &model_conf);
        let psi_star = solve_poisson(&mut psi, &source, &model_conf, &bc);
        let diff = &psi_star + &(-1.0 * &psi);
        psi = &psi + &(model_conf.relax * diff);
        psi.fill_bc(&bc);
        progress_bar.increment(1);
    }
    progress_bar.finish();

    io::write_vtk(&"./examples/magnetostatics/", &"ferro_droplet", model_conf.n_iter, &psi, &output, &progress_bar);
    println!("{}", "Done.".green().bold());

}