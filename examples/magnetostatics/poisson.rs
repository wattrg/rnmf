
use super::mesh::cartesian2d::*;
use super::boundary_conditions::{BcType, BCs, ComponentBCs};
use crate::config::Config;
use super::model::UserModel;
use rnmf::{Real, RealVec2};

fn get_laplace(phi: &CartesianDataFrame2D<Real>, model: &Config<UserModel>) -> CartesianDataFrame2D<Real> {
    let bc = BCs::new(
        vec![
            ComponentBCs::new(
                vec![BcType::Dirichlet(0.0), BcType::Dirichlet(0.0)],
                vec![BcType::Dirichlet(0.0), BcType::Dirichlet(0.0)]
            )
        ]
    );
    let mut lap = CartesianDataFrame2D::new_from(&phi.underlying_mesh, bc, 1, 1);
    let dx = &phi.underlying_mesh.dx;

    for ((i,j,_), lap) in lap.iter_mut().enumerate_index(){
        *lap = phi[(i+1,j  ,0)]*(mu(i+1,j  ,dx,model) + mu(i,j,dx,model))/2.0
                        + phi[(i-1,j  ,0)]*(mu(i-1,j  ,dx,model) + mu(i,j,dx,model))/2.0
                        + phi[(i  ,j+1,0)]*(mu(i  ,j+1,dx,model) + mu(i,j,dx,model))/2.0
                        + phi[(i  ,j-1,0)]*(mu(i  ,j-1,dx,model) + mu(i,j,dx,model))/2.0
                        - (phi[(i  ,j  ,0)]*(mu(i+1,j  ,dx,model) + mu(i-1,j,dx,model)
                                +mu(i,j+1,dx,model)+mu(i,j-1,dx,model)+4.0*mu(i,j,dx,model)))/2.0;
    }
    
    lap.fill_bc();
    lap
}

pub fn get_source(phi: &CartesianDataFrame2D<Real>, config: &Config<UserModel>) -> CartesianDataFrame2D<Real> {
    let mut source = CartesianDataFrame2D::new_from(&phi.underlying_mesh, phi.bc.clone(), 1, 1);
    for ((i,j,_), src) in &mut source.iter_mut().enumerate_index(){
        *src = -(1.0 - 1.0/config.model.relax)/mu(i,j,&phi.underlying_mesh.dx, config);
    }
    &source * &get_laplace(phi, config)
}


fn is_inside(i: isize, j: isize, dx: &RealVec2, config: &Config<UserModel>) -> bool{
    let x_pos = (i as Real + 0.5) * dx[0];
    let y_pos = (j as Real + 0.5) * dx[1];

    let x_dist_squared = (x_pos - config.model.bubble_centre[0]).powf(2.0);
    let y_dist_squared = (y_pos - config.model.bubble_centre[1]).powf(2.0);

    x_dist_squared + y_dist_squared <= config.model.bubble_radius*config.model.bubble_radius
}

pub fn mu(i: isize, j: isize, dx: &RealVec2, config: &Config<UserModel>) -> Real {
    if is_inside(i, j, dx, config){
        config.model.mu[0]
    }
    else {
        config.model.mu[1]
    }
}


pub fn solve_poisson(psi_init: &CartesianDataFrame2D<Real>, 
                     rhs: &CartesianDataFrame2D<Real>, 
                     config: &Config<UserModel>) -> CartesianDataFrame2D<Real>{

    let mut psi = psi_init.clone();
    let dx2 = psi.underlying_mesh.dx[0].powf(2.0);
    let dy2 = psi.underlying_mesh.dx[1].powf(2.0);
    let alpha = dx2 * dy2 / (2.0 * (dx2 + dy2));
    let beta = dy2 / (2.0 * (dx2 + dy2));
    let gamma = dx2 / (2.0 * (dx2 + dy2));

    // since we're directly updating psi with itself, we can't use the IterMut trait
    // for the moment we just have to use regular loops
    for _ in 0..config.model.n_sub_iter {
        for i in 0..psi.underlying_mesh.n[0] as isize{
            for j in 0..psi.underlying_mesh.n[1] as isize {
                psi[(i,j,0)] = alpha * rhs[(i,j,0)] + 
                                 beta * (psi[(i+1,j,0)] + psi[(i-1,j,0)]) + 
                                 gamma * (psi[(i,j+1,0)] + psi[(i,j-1,0)]);
            }
        }
        psi.fill_bc();
    }
    psi

}

pub fn sum(psi: &CartesianDataFrame2D<Real>)->Real{
    let mut psi_sum: Real = 0.0;
    for val in psi.iter(){
        psi_sum += val.abs();
    }
    psi_sum
}

