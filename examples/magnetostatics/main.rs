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
    let commit = &String::from_utf8(
        Command::new("git").args(&["rev-parse", "HEAD"])
                           .output()
                           .unwrap()
                           .stdout).unwrap()[0..6];

    println!("{}", format!("rnmf-{}: version {}", commit, version).green().bold());
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


    let model_instance = model::UserModel::new();
    let conf = config::read_lua(&args[1], model_instance).unwrap();



    // set up hash map for outputting variables
    let mut output: HashMap<String, io::OutputCallBack> = HashMap::new();
    output.insert("psi".to_string(), model::get_psi);
    output.insert("h".to_string(), model::get_h);

    // create the mesh 
    let mesh = CartesianMesh2D::new(
        [0.0, 0.0],                                   // lo corner
        [conf.geom.length[0], conf.geom.length[1]],   // hi corner
        [conf.geom.n_cells[0], conf.geom.n_cells[1]]  // number of cells
    );


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

    io::write_vtk(&"./examples/magnetostatics/", &"ferro_droplet", 0, &psi, &output, &progress_bar);

    // begin solving
    for _ in 0..conf.model.n_iter{
        let source = get_source(&psi, &conf);
        let psi_star = solve_poisson(&mut psi, &source, &conf, &bc);
        let diff = &psi_star + &(-1.0 * &psi);
        psi = &psi + &(conf.model.relax * diff);
        psi.fill_bc(&bc);
        progress_bar.increment(1);
    }
    progress_bar.finish();

    io::write_vtk(&"./examples/magnetostatics/", &"ferro_droplet", conf.model.n_iter, &psi, &output, &progress_bar);
    println!("{}", "Done.".green().bold());

}