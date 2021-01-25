use rlua::{UserData, Context};
use crate::{Real, RealVec2};
use std::convert::{TryFrom};
use crate::mesh::cartesian2d::{CartesianDataFrame2D,IndexEnumerable};
use crate::config::UserConfig;
use crate::io::{VtkData};

#[derive(Clone, Debug)]
pub struct UserModel{
    pub h_far: RealVec2,
    pub relax: Real,
    pub n_iter: usize,
    pub n_sub_iter: usize,
    pub bubble_centre: RealVec2,
    pub bubble_radius: Real,
    pub mu: RealVec2,
}
impl UserModel{
    pub fn new()->Self{
        Self{
            h_far: RealVec2([0.0, 0.0]),
            bubble_centre: RealVec2([0.0, 0.0]),
            mu: RealVec2([0.0, 0.0]),
            n_iter: 0,
            n_sub_iter: 0,
            relax: 0.0,
            bubble_radius: 0.0,
        }
    }
}
impl UserData for UserModel {}

impl UserConfig for UserModel{
    fn lua_constructor(self, lua_ctx: Context)->rlua::Function{
        lua_ctx.create_function(|_,(h_far,bubble_centre,mu,n_iter,n_sub_iter,relax,bubble_radius): (RealVec2,RealVec2,RealVec2,usize,usize,Real,Real)|
            Ok(UserModel{
                h_far,
                bubble_centre,
                mu,
                n_iter,
                n_sub_iter,
                relax,
                bubble_radius,
            })).expect("failed creating user model lua constructor")
    }
}

pub fn get_psi(data: &CartesianDataFrame2D) -> Result<VtkData, String> {
    Ok(VtkData::try_from(data.clone())?)
}

pub fn get_h(data: &CartesianDataFrame2D) -> Result<VtkData, String> {
    let mut h_field = CartesianDataFrame2D::new_from(&data.underlying_mesh, 2, data.n_ghost);
    for ((i,j,n), h) in h_field.into_iter().enumerate_index(){
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
