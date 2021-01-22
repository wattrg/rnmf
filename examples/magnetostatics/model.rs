use rlua::{Lua};
use std::fs;
use crate::{Real, RealVec2, UIntVec2};
use std::convert::{TryInto};

#[derive(Clone)]
pub struct ModelConfig {
    pub h_far: RealVec2,
    pub length: RealVec2,
    pub n_cells: UIntVec2,
    pub relax: Real,
    pub n_iter: usize,
    pub n_sub_iter: usize,
    pub bubble_centre: RealVec2,
    pub bubble_radius: Real,
    pub mu: RealVec2,
}


// this should be done elsewhere
pub fn read_lua(lua_loc: &str) -> Result<ModelConfig, std::io::Error>{
    // load the contents of the file
    let lua_file = fs::read_to_string(lua_loc).expect("Could not read lua file");

    // create the lua struct
    let lua = Lua::new();

    // initialise the configuration struct
    let mut model_conf = ModelConfig{
        h_far: [0.0, 1.0],
        length: [1.0, 1.0],
        n_cells: [10, 10],
        relax: 0.4,
        n_iter: 0,
        bubble_centre: [0.5,0.5],
        bubble_radius: 0.1,
        mu: [1.0, 1.0],
        n_sub_iter: 0,
    };

    // create lua context to interface with lua
    lua.context(|lua_ctx|{
        // access lua globals
        let globals = lua_ctx.globals();

        // execute the lua script
        lua_ctx.load(&lua_file)
               .exec()
               .expect("failed executing lua file");

        // set configuration data for rust
        model_conf.h_far = globals.get::<_,Vec<Real>>("H_far")
                                .expect("failed setting dimensionality")[0..2]
                                .try_into()
                                .expect("wrong number of elements in H_far");
        println!("Set h_far to {:?}", model_conf.h_far);

        model_conf.length = globals.get::<_,Vec<Real>>("length")
                                .expect("failed setting domain length")[0..2]
                                .try_into()
                                .expect("wrong number of elements in 'length'");
        println!("Set length to {:?}", model_conf.length);

        model_conf.n_cells = globals.get::<_,Vec<usize>>("n_cells")
                                .expect("failed setting number of cells")[0..2]
                                .try_into()
                                .expect("wrong number of elements in 'n_cells'");
        println!("Set n_cells to {:?}", model_conf.n_cells);

        model_conf.relax = globals.get::<_,Real>("relax")
                                .expect("failed setting relaxation coefficient");
        println!("Set relaxation coefficient to {}", model_conf.relax);

        model_conf.n_iter = globals.get::<_,usize>("n_iter")
                                .expect("failed setting number of iterations");
        println!("Set n_iter to {}", model_conf.n_iter);

        model_conf.n_sub_iter = globals.get::<_,usize>("n_sub_iter")
                                    .expect("failed setting n_sub_iter");
        println!("Set n_sub_iter to {}", model_conf.n_sub_iter);

        model_conf.bubble_centre = globals.get::<_,Vec<Real>>("bubble_centre")
                                        .expect("failed setting bubble_centre")[0..2]
                                        .try_into()
                                        .expect("wrong number of elements in 'bubble_centre'");
        println!("Set bubble centre to {:?}", model_conf.bubble_centre);

        model_conf.bubble_radius = globals.get::<_,Real>("bubble_radius")
                                        .expect("failed setting bubble radius");
        println!("Set bubble radius to {}", model_conf.bubble_radius);

        model_conf.mu = globals.get::<_,Vec<Real>>("mu")
                                .expect("failed setting mu")[0..2]
                                .try_into()
                                .expect("wrong number of elements in 'mu'");
        println!("Set mu_1 to {:?}", model_conf.mu);
        println!();
    });
    Ok(model_conf)
}