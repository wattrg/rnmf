use rlua::{UserData, Context};
use crate::{Real, RealVec2};
use std::convert::{TryFrom};
use crate::mesh::cartesian2d::{CartesianDataFrame2D,IndexEnumerable};
use crate::config::UserConfig;
use crate::io::{VtkData};
use std::collections::HashMap;

#[derive(Clone, Debug)]
pub struct UserModel{
    pub h_far: RealVec2,
    pub relax: Real,
    pub n_iter: usize,
    pub n_sub_iter: usize,
    pub bubble_centre: RealVec2,
    pub bubble_radius: Real,
    pub mu: RealVec2,
    pub tol: Real,
}
impl UserData for UserModel {}

impl UserConfig for UserModel{
    fn new()->Self{
        Self{
            h_far: RealVec2([0.0, 0.0]),
            bubble_centre: RealVec2([0.0, 0.0]),
            mu: RealVec2([0.0, 0.0]),
            n_iter: 0,
            n_sub_iter: 0,
            relax: 0.0,
            bubble_radius: 0.0,
            tol: 0.4,
        }
    }

    fn lua_constructor(lua_ctx: Context)->rlua::Function{
        lua_ctx.create_function(|_,model: rlua::Table|
            Ok(UserModel{
                h_far: model.get::<_,RealVec2>("H_far")
                            .expect("failed reading 'H_far'"),
                bubble_centre: model.get::<_,RealVec2>("bubble_centre")
                                    .expect("failed reading 'bubble centre'"),
                mu: model.get::<_,RealVec2>("mu")
                         .expect("failed reading 'mu'"),
                n_iter: model.get::<_,usize>("n_iter")
                             .expect("failed reading 'n_iter'"),
                n_sub_iter: model.get::<_,usize>("n_sub_iter")
                                 .expect("failed reading 'n_sub_iter'"),
                relax: model.get::<_,Real>("relax")
                            .expect("failed reading 'relax'"),
                bubble_radius: model.get::<_,Real>("bubble_radius")
                                    .expect("failed reading 'bubble_radius'"),
                tol: model.get::<_,Real>("tol")
                          .expect("failed reading 'tol'")
            })
        ).expect("failed creating user model lua ")
    }
}

pub fn output_callbacks()-> crate::io::OutputCallBackHashMap<Real>{
    // set up hash map for outputting variables
    let mut output: crate::io::OutputCallBackHashMap<Real> = HashMap::new();
    output.insert("psi".to_string(), get_psi);
    output.insert("h".to_string(), get_h);
    output
}

pub fn get_psi(data: &CartesianDataFrame2D<Real>) -> Result<VtkData<f64>, String> {
    Ok(VtkData::try_from(data.clone())?)
}

pub fn get_h(data: &CartesianDataFrame2D<Real>) -> Result<VtkData<Real>, String> {
    let mut h_field = CartesianDataFrame2D::new_from(&data.underlying_mesh, data.bc.clone(), 2, data.n_ghost);
    for ((i,j,n), h) in &mut h_field.iter_mut().enumerate_index(){
        match n {
            0 => {
                *h = -(data[(i+1,j,0)] - data[(i-1,j,0)])/(2.0*data.underlying_mesh.dx[0]);
            }
            1 => {
                *h = -(data[(i,j+1,0)] - data[(i,j-1,0)])/(2.0*data.underlying_mesh.dx[1]);
            }
            _ => {
                return Err(format!("Component {} out of bounds of data with {} components", n, h_field.n_comp));
            }
        }
    }
    Ok(VtkData::try_from(h_field)?)
}
