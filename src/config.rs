use rlua::{Lua, UserData};
use std::fs;
use crate::*;
use crate::mesh::cartesian2d::CartesianDataFrame2D;
use colored::*;

pub trait UserConfig {}

pub struct Actions<T: UserConfig>{
    pub actions: Vec<Action<T>>
}
impl <T: UserConfig> Actions<T>{
    pub fn new()->Actions<T>{
        Actions{
            actions: Vec::new(),
        }
    }
}

pub struct Action<T: UserConfig>{
    pub name: String,
    pub pre_action: Option<fn(&Config<T>)>,
    pub action: Option<fn(&CartesianDataFrame2D, &Config<T>)>,
    pub iters: usize,
    pub stop: Option<f64>,
}
impl <T: UserConfig> Action<T>{
    pub fn new()->Action<T>{
        Action{
            name: String::from(""),
            pre_action: Option::None,
            action: Option::None,
            iters: 0,
            stop: Option::None,
        }
    }
}

pub struct Config<T: UserConfig>{
    pub geom: GeomConfig,
    pub model: T,
    pub actions: Actions<T>,
}
impl <T: UserConfig> Config<T>{
    pub fn new(user_config: T)->Config<T>{
        Config{
            geom: GeomConfig::new(),
            model: user_config,
            actions: Actions::new(),
        }
    }
}

pub struct GeomConfig{
    pub dim: usize,
    pub length: RealVec2,
    pub n_cells: UIntVec2,
}
impl GeomConfig{
    pub fn new() -> GeomConfig{
        GeomConfig{
            dim: 0,
            length: RealVec2([0.0, 0.0]),
            n_cells: UIntVec2([0, 0]),
        }
    }
}

impl UserData for &RnmfType {}
impl UserData for RnmfType {}
impl UserData for &mut &RnmfType  {}
impl UserData for RealVec2 {}
impl UserData for UIntVec2 {}
impl UserData for IntVec2 {}
impl UserData for RealVec3 {}
impl UserData for UIntVec3 {}
impl UserData for IntVec3 {}

/// function which executes a lua file (located at lua_loc), and returns the configuration
pub fn read_lua<T>(lua_loc: &str, user_model: T) -> Result<Config<T>, std::io::Error>
    where 
        T: UserConfig,
{
    let mut conf = Config::new(user_model);

    // load the contents of the file
    let lua_file = fs::read_to_string(lua_loc).expect("Could not read lua file");

    // create the lua struct
    let lua = Lua::new();


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

        // read geometry data from the lua file
        conf.geom.dim = globals.get::<_,usize>("dim")
                          .expect("failed setting dimensionality");

        match conf.geom.dim{
            2 => {
                conf.geom.length = globals.get::<_,RealVec2>("length")
                                          .unwrap();
                conf.geom.n_cells = globals.get::<_,UIntVec2>("n_cells")
                                           .unwrap();
            }
            _ => {
                println!("{} {} not yet supported. Exiting", "Error:".red().bold(), conf.geom.dim);
            }
        }

    });
    Ok(conf)
}
