mod poisson;

use rnmf::*;
use boundary_conditions::*;
use mesh::cartesian::*;
use std::env;
use poisson::*;
use io::{RnmfProgressBar};

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

    // some config data which will eventually go into the lua file
    let h_far = vec![0.0, 1.0];
    let l_x = 1.0;
    let l_y = 1.0;
    let nx  = 5;
    let ny  = 5;
    let epsilon = 0.005;

    // create the mesh 
    let mesh = CartesianMesh::new(vec![0.0, 0.0], vec![l_x, l_y], vec![nx, ny], conf.dim);

    let bc = BCs::new(vec![
        ComponentBCs::new(
            vec![BCType::Neumann(h_far[0]/mesh.dx[0]), BCType::Neumann(-h_far[1]/mesh.dx[1])], // lo bc
            vec![BCType::Neumann(h_far[0]/mesh.dx[0]), BCType::Neumann(-h_far[1]/mesh.dx[1])]  // hi bc
        )
    ]);

    // create the data frame
    let mut psi = CartesianDataFrame::new_from(&mesh, 1, 1);

    // fill the initial condition
    psi.fill_ic(|x,y,_,_| -h_far[0] * x - h_far[1] * y);

    // fill the boundary conditions
    psi.fill_bc(bc);

    // begin solving

    println!("Overall progress:");
    let mut progress_bar = io::ProgressBar::create(1000);

    let h_field = -1.0*poisson::gradient(&psi);
    let mut b_i = mul_mu(&h_field);

    //let mut iter = 1;

    for _ in 0..1000{
        let lap_psi_new = divide_mu(&((1.0 - 1.0/epsilon) * divergence(&b_i)));
        solve_poisson(&mut psi, &lap_psi_new);
        let b_star = mul_mu(&(-1.0*poisson::gradient(&psi)));
        let diff = &b_star + &(-1.0 * &b_i);
        b_i = &b_i + &(epsilon * diff);
        let bc = BCs::new(
            vec![
                // x component
                ComponentBCs::new(
                    vec![BCType::Dirichlet(0.0), BCType::Dirichlet(0.0)], // lo bc
                    vec![BCType::Dirichlet(0.0), BCType::Dirichlet(0.0)]  // hi bc
                ),
                // y component
                ComponentBCs::new(
                    vec![BCType::Dirichlet(1.0), BCType::Dirichlet(1.0)], // lo bc
                    vec![BCType::Dirichlet(1.0), BCType::Dirichlet(1.0)]  // hi bc
                )
            ],
        );
        b_i.fill_bc(bc);
        // if iter % 10 == 0 {
        // }
       // iter += 1;
        progress_bar.increment(1);
    }
    progress_bar.finish();

    psi = io::write_csv(&"./examples/magnetostatics/test", psi);
    
    println!("Done.");

}