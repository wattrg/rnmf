
use super::mesh::cartesian::*;
use super::boundary_conditions::*;
use super::model::ModelConfig;
use rnmf::Real;

fn get_laplace(phi: &CartesianDataFrame, model: &ModelConfig) -> CartesianDataFrame {
    let mut lap = CartesianDataFrame::new_from(&phi.underlying_mesh, 1, 1);
    let dx = &phi.underlying_mesh.dx;
    for i in 0..phi.underlying_mesh.n[0] as isize {
        for j in 0..phi.underlying_mesh.n[1] as isize {
            lap[(i,j,0,0)] = phi[(i+1,j  ,0,0)]*(mu(i+1,j  ,dx,model) + mu(i,j,dx,model))/2.0
                           + phi[(i-1,j  ,0,0)]*(mu(i-1,j  ,dx,model) + mu(i,j,dx,model))/2.0
                           + phi[(i  ,j+1,0,0)]*(mu(i  ,j+1,dx,model) + mu(i,j,dx,model))/2.0
                           + phi[(i  ,j-1,0,0)]*(mu(i  ,j-1,dx,model) + mu(i,j,dx,model))/2.0
                           - phi[(i  ,j  ,0,0)]*(mu(i+1,j  ,dx,model) + mu(i-1,j,dx,model)+mu(i,j+1,dx,model)+mu(i,j-1,dx,model)+4.0*mu(i,j,dx,model))/2.0;
        }
    }
    let bc = BCs::new(
        vec![
            ComponentBCs::new(
                vec![BCType::Dirichlet(0.0), BCType::Dirichlet(0.0)],
                vec![BCType::Dirichlet(0.0), BCType::Dirichlet(0.0)]
            )
        ]
    );
    lap.fill_bc(&bc);
    lap
}

pub fn get_source(phi: &CartesianDataFrame, model: &ModelConfig) -> CartesianDataFrame {
    let mut source = CartesianDataFrame::new_from(&phi.underlying_mesh, 1, 1);
    for ((i,j,_,_), src) in source.into_iter().enumerate_index(){
        *src = -(1.0 - 1.0/model.relax)/mu(i,j,&phi.underlying_mesh.dx, model);
    }
    &source * &get_laplace(phi, model)
}


fn is_inside(i: isize, j: isize, dx: &[Real], model: &ModelConfig) -> bool{
    let x_pos = (i as Real + 0.5) * dx[0];
    let y_pos = (j as Real + 0.5) * dx[1];

    let x_dist = (x_pos - model.bubble_centre[0]).abs();
    let y_dist = (y_pos - model.bubble_centre[1]).abs();

    x_dist * x_dist + y_dist * y_dist <= model.bubble_radius*model.bubble_radius
}

fn mu(i: isize, j: isize, dx: &[Real], model: & ModelConfig) -> Real {
    if is_inside(i, j, dx, model){
        model.mu[0]
    }
    else {
        model.mu[1]
    }
}


pub fn solve_poisson(psi_init: &CartesianDataFrame, 
                     rhs: &CartesianDataFrame, 
                     model: &ModelConfig,
                     bc: &BCs) -> CartesianDataFrame{

    let mut psi = psi_init.clone();
    let dx2 = psi.underlying_mesh.dx[0].powf(2.0);
    let dy2 = psi.underlying_mesh.dx[1].powf(2.0);
    let alpha = dx2 * dy2 / (2.0 * (dx2 + dy2));
    let beta = dy2 / (2.0 * (dx2 + dy2));
    let gamma = dx2 / (2.0 * (dx2 + dy2));


    for _ in 0..model.n_sub_iter {
        for i in 0..psi.underlying_mesh.n[0] as isize{
            for j in 0..psi.underlying_mesh.n[1] as isize {
                psi[(i,j,0,0)] = alpha * rhs[(i,j,0,0)] + 
                                 beta * (psi[(i+1,j,0,0)] + psi[(i-1,j,0,0)]) + 
                                 gamma * (psi[(i,j+1,0,0)] + psi[(i,j-1,0,0)]);
            }
        }
        psi.fill_bc(&bc);
    }
    psi

}

