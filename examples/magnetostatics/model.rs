use rlua::{Lua};
use std::fs;
use crate::{Real, RealVec2, UIntVec2, IntVec2};
use std::convert::{TryFrom};
use crate::mesh::cartesian2d::*;
use crate::io::*;
//use crate::config::*;


pub struct ModelConfig {
    pub h_far: RealVec2,
    pub relax: Real,
    pub n_iter: usize,
    pub n_sub_iter: usize,
    pub bubble_centre: RealVec2,
    pub bubble_radius: Real,
    pub mu: RealVec2,
}


/// struct which stores all the geometric configuration for the simulation
pub struct GeomConfig {
    pub dim: usize,
    pub length: RealVec2,
    pub n_cells: UIntVec2,
}
pub struct Config {
    pub geom: GeomConfig,
    pub model: ModelConfig,
    //pub actions: Option<Actions<f64>>,
}



// this should be done elsewhere
pub fn read_lua(lua_loc: &str) -> Result<Config, std::io::Error>{
    // load the contents of the file
    let lua_file = fs::read_to_string(lua_loc).expect("Could not read lua file");

    // create the lua struct
    let lua = Lua::new();

    // initialise the configuration struct
    let mut model_conf = ModelConfig{
        h_far: RealVec2([0.0, 1.0]),
        relax: 0.4,
        n_iter: 0,
        bubble_centre: RealVec2([0.5,0.5]),
        bubble_radius: 0.1,
        mu: RealVec2([1.0, 1.0]),
        n_sub_iter: 0,
    };

    let mut geom_conf = GeomConfig{
        length: RealVec2([1.0, 1.0]),
        n_cells: UIntVec2([10, 10]),
        dim: 0,
    };

    // create lua context to interface with lua
    lua.context(|lua_ctx|{
        // access lua globals
        let globals = lua_ctx.globals();

        // create constructors for data types
        let realvec2_constructor = lua_ctx.create_function(|_,(x,y): (Real, Real)| Ok(RealVec2([x, y])))
                                          .expect("Failed creating 'RealVec2' object from lua");
        let uintvec2_constructor = lua_ctx.create_function(|_,(x,y): (usize, usize)| Ok(UIntVec2([x, y])))
                                          .expect("Failed creating 'UIntVec' object from lua");
        let intvec2_constructor = lua_ctx.create_function(|_,(x,y): (isize, isize)| Ok(IntVec2([x, y])))
                                          .expect("Failed creating 'UIntVec' object from lua");

        globals.set("RealVec2", realvec2_constructor).unwrap();
        globals.set("UIntVec2", uintvec2_constructor).unwrap();
        globals.set("IntVec2", intvec2_constructor).unwrap();

        // execute the lua script
        lua_ctx.load(&lua_file)
               .exec()
               .expect("failed executing lua file");

        // set configuration data for rust
        model_conf.h_far = globals.get::<_,RealVec2>("H_far")
                                .expect("failed setting H_far");
        println!("Set h_far to {:?}", model_conf.h_far);

        geom_conf.length = globals.get::<_,RealVec2>("length")
                                .expect("failed setting domain length");
        println!("Set length to {:?}", geom_conf.length);

        geom_conf.n_cells = globals.get::<_,UIntVec2>("n_cells")
                                .expect("failed setting number of cells");
        println!("Set n_cells to {:?}", geom_conf.n_cells);

        model_conf.relax = globals.get::<_,Real>("relax")
                                .expect("failed setting relaxation coefficient");
        println!("Set relaxation coefficient to {}", model_conf.relax);

        model_conf.n_iter = globals.get::<_,usize>("n_iter")
                                .expect("failed setting number of iterations");
        println!("Set n_iter to {}", model_conf.n_iter);

        model_conf.n_sub_iter = globals.get::<_,usize>("n_sub_iter")
                                    .expect("failed setting n_sub_iter");
        println!("Set n_sub_iter to {}", model_conf.n_sub_iter);

        model_conf.bubble_centre = globals.get::<_,RealVec2>("bubble_centre")
                                        .expect("failed setting bubble_centre");
        println!("Set bubble centre to {:?}", model_conf.bubble_centre);

        model_conf.bubble_radius = globals.get::<_,Real>("bubble_radius")
                                        .expect("failed setting bubble radius");
        println!("Set bubble radius to {}", model_conf.bubble_radius);

        model_conf.mu = globals.get::<_,RealVec2>("mu")
                                .expect("failed setting mu");
        println!("Set mu_1 to {:?}", model_conf.mu);
        println!();
    });
    Ok(Config{
        model: model_conf,
        geom: geom_conf,
        //actions: Option::None,
    })
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
