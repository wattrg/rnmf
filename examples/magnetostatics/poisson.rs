
use super::mesh::cartesian::*;
use super::boundary_conditions::*;


pub fn gradient(phi: &CartesianDataFrame) -> CartesianDataFrame {
    let mut grad = CartesianDataFrame::new_from(&phi.underlying_mesh, 2, 1);
    for i in 0..phi.underlying_mesh.n[0] as isize{
        for j in 0..phi.underlying_mesh.n[1] as isize{
            grad[(i,j,0,0)] = (phi[(i+1,j,0,0)] - phi[(i-1,j,0,0)])/(2.0*phi.underlying_mesh.dx[0]);
            grad[(i,j,0,1)] = (phi[(i,j+1,0,0)] - phi[(i,j-1,0,0)])/(2.0*phi.underlying_mesh.dx[1]);
        }
    }
    let bc = BCs::new(
        vec![
            // x component
            ComponentBCs::new(
                vec![BCType::Dirichlet(0.0), BCType::Dirichlet(0.0)], // lo bc
                vec![BCType::Dirichlet(0.0), BCType::Dirichlet(0.0)]  // hi bc
            ),
            // y component
            ComponentBCs::new(
                vec![BCType::Dirichlet(-1.0), BCType::Dirichlet(-1.0)], // lo bc
                vec![BCType::Dirichlet(-1.0), BCType::Dirichlet(-1.0)]  // hi bc
            )
        ],
    );
    grad.fill_bc(bc);
    grad
}

// pub fn laplace(phi: &CartesianDataFrame) -> CartesianDataFrame {
//     let mut lap = CartesianDataFrame::new_from(&phi.underlying_mesh, 1, 0);
//     for i in 0..phi.underlying_mesh.n[0] as isize {
//         for j in 0..phi.underlying_mesh.n[1] as isize {
//             //lap[(i,j,0,0)] = 
//         }
//     }
//     lap
// }

pub fn divergence(vec: &CartesianDataFrame) -> CartesianDataFrame {
    let mut div = CartesianDataFrame::new_from(&vec.underlying_mesh, 1, 1);
    for i in 0..vec.underlying_mesh.n[0] as isize{
        for j in 0..vec.underlying_mesh.n[1] as isize {
            let d_x = (vec[(i+1,j,0,0)] - vec[(i-1,j,0,0)])/(2.0*vec.underlying_mesh.dx[0]);
            let d_y = (vec[(i,j+1,0,1)] - vec[(i,j-1,0,1)])/(2.0*vec.underlying_mesh.dx[1]);
            div[(i,j,0,0)] = d_x + d_y;
        }
    }
    div
}

fn is_inside(i: isize, j: isize, dx: &[f64], r: f64, cntr: &[f64]) -> bool{
    let x_pos = (i as f64 + 0.5) * dx[0];
    let y_pos = (j as f64 + 0.5) * dx[1];

    let x_dist = (x_pos - cntr[0]).abs();
    let y_dist = (y_pos - cntr[1]).abs();

    x_dist * x_dist + y_dist * y_dist <= r*r
}

pub fn mu(i: isize, j: isize, dx: &[f64], r: f64, cntr: &[f64]) -> f64 {
    if is_inside(i, j, dx, r, cntr){
        3.0
    }
    else {
        1.0
    }
}

pub fn mul_mu (h_field: &CartesianDataFrame) -> CartesianDataFrame {
    let mut b_field = CartesianDataFrame::new_from(&h_field.underlying_mesh, h_field.n_comp, 1);
    for i in 0..h_field.underlying_mesh.n[0] as isize{
        for j in 0..h_field.underlying_mesh.n[1] as isize {
            for n in 0..h_field.n_comp {
                b_field[(i,j,0,n)] = h_field[(i,j,0,n)] * mu(i, j, 
                                                             &b_field.underlying_mesh.dx,
                                                             0.5, 
                                                             &b_field.underlying_mesh.hi
                                     );
            }
        }
    }
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
    b_field.fill_bc(bc);
    b_field
}

pub fn divide_mu(scalar: &CartesianDataFrame) -> CartesianDataFrame {
    let mut ret = CartesianDataFrame::new_from(&scalar.underlying_mesh, scalar.n_comp, 1);
    
    for i in 0..scalar.underlying_mesh.n[0] as isize{
        for j in 0..scalar.underlying_mesh.n[1] as isize {
            for n in 0..scalar.n_comp {
                ret[(i,j,0,n)] = scalar[(i,j,0,n)] /  mu(i, j, 
                                                             &scalar.underlying_mesh.dx,
                                                             0.5, 
                                                             &scalar.underlying_mesh.hi
                                     );
            }
        }
    }
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
    ret.fill_bc(bc);
    ret

}


pub fn solve_poisson(psi: &mut CartesianDataFrame, rhs: &CartesianDataFrame){
    let dx2 = psi.underlying_mesh.dx[0].powf(2.0);
    let dy2 = psi.underlying_mesh.dx[1].powf(2.0);
    let alpha = dx2 * dy2 / (2.0 * (dx2 + dy2));
    let beta = dy2 / (2.0 * (dx2 + dy2));
    let gamma = dx2 / (2.0 * (dx2 + dy2));


    for _ in 0..1000 {
        for i in 0..psi.underlying_mesh.n[0] as isize{
            for j in 0..psi.underlying_mesh.n[1] as isize {
                psi[(i,j,0,0)] = alpha * rhs[(i,j,0,0)] + beta * (psi[(i+1,j,0,0)] + psi[(i+1,j,0,0)]) + gamma * (psi[(i,j+1,0,0)] + psi[(i,j-1,0,0)]);
            }
        }
        let bc = BCs::new(vec![
            ComponentBCs::new(
                vec![BCType::Neumann(0.0), BCType::Neumann(-1.0/psi.underlying_mesh.dx[1])], // lo bc
                vec![BCType::Neumann(0.0), BCType::Neumann(-1.0/psi.underlying_mesh.dx[1])]  // hi bc
            )
        ]);
        psi.fill_bc(bc);
    }

}

