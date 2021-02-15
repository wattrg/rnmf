use rlua::{Lua, UserData, Context};
use std::fs;
use crate::*;
use colored::*;
use std::process::Command;
use std::collections::HashMap;

pub trait UserConfig: UserData + Clone {
    fn new()->Self;
    fn lua_constructor(lua_ctx: Context)->rlua::Function;
}

#[derive(Clone)]
pub struct Actions{
    pub actions: Vec<Action>
}
impl UserData for Actions{}
impl UserConfig for Actions{
    fn new()->Actions{
        Actions{
            actions: Vec::new(),
        }
    }

    fn lua_constructor(lua_ctx: Context)->rlua::Function{
        lua_ctx.create_function(|_,actions: Action|
            Ok(
                //Actions {actions: actions.get::<_,Vec<Action>>(1).expect("failed reading actions"),}
                Actions {actions: vec![actions]}
            )
        ).expect("failed reading Actions from lua file")
    }

}

#[derive(Clone)]
pub struct Action{
    pub name: String,
    //pub pre_action: Option<fn(&Config<T>)>,
    //pub action: Option<fn(&CartesianDataFrame2D, &Config<T>)>,
    pub action: String,
    pub iters: usize,
    //pub stop: Option<f64>,
}
impl UserConfig for Action {
    fn new()->Self{
        Self{
            name: String::from(""),
            action: String::from(""),
            iters: 0,
        }
    }

    fn lua_constructor(lua_ctx: Context)->rlua::Function{
        lua_ctx.create_function(|_,action: rlua::Table|
            Ok(
                Action{
                    name: action.get::<_,String>("name").expect("failed reading action name"),
                    action: action.get::<_,String>("action").expect("failed reading action"),
                    iters: action.get::<_,usize>("iterations")
                                 .expect("failed reading number of iterations for action"),
                })
        ).expect("failed creating action from lua file ")
    }
}

impl  UserData for Action{}



pub struct Config<T: UserConfig>{
    pub geom: GeomConfig,
    pub model: T,
    pub actions: Actions,
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

impl UserData for RealVec2 {}
impl UserData for UIntVec2 {}
impl UserData for IntVec2 {}
impl UserData for RealVec3 {}
impl UserData for UIntVec3 {}
impl UserData for IntVec3 {}

pub fn init<T: 'static>(args: Vec<String>, user_model: T, out_cb:fn()->HashMap<String, crate::io::OutputCallBack>)
    ->Result<(Config<T>, HashMap<String,crate::io::OutputCallBack>), std::io::Error>
    where 
        T: UserConfig    
{
    let version = env!("CARGO_PKG_VERSION");
    let commit = &String::from_utf8(
        Command::new("git").args(&["rev-parse", "HEAD"])
                           .output()
                           .unwrap()
                           .stdout).unwrap()[0..7];

    println!("{}", format!("rnmf {}-{}", version, commit).green().bold());
    if cfg!(feature = "disable_double") {
        println!("using single precision");
    }
    else{
        println!("using double precision");
    }
    println!("Hello from the magnetistatics example!");

    if args.len() < 1 {
        println!("{} Location of lua configuration script not given.", "Error:".red());
        panic!();
    }

    Ok((read_lua(&args[1], user_model)?, out_cb()))


}

/// function which executes a lua file (located at lua_loc), and returns the configuration
pub fn read_lua<T: 'static>(lua_loc: &str, user_model: T) -> Result<Config<T>, std::io::Error>
    where 
        T: UserConfig,
{
    let mut conf = Config::new(user_model.clone());

    // load the contents of the file
    let lua_file = fs::read_to_string(lua_loc).expect("Could not read lua file");

    // create the lua struct
    let lua = Lua::new();


    // create lua context to interface with lua
    lua.context(|lua_ctx|{
        // access lua globals
        let globals = lua_ctx.globals();

        // create constructors for data types
        let realvec2_constructor = lua_ctx.create_function(|_,(x,y): (Real, Real)| 
            Ok(RealVec2([x, y]))
        ).expect("Failed creating 'RealVec2' constructor for lua");
        let uintvec2_constructor = lua_ctx.create_function(|_,(x,y): (usize, usize)| 
            Ok(UIntVec2([x, y]))
        ).expect("Failed creating 'UIntVec' constructor for lua");
        let intvec2_constructor = lua_ctx.create_function(|_,(x,y): (isize, isize)| 
            Ok(IntVec2([x, y]))
        ).expect("Failed creating 'UIntVec' constructor for lua");

        let user_model_constructor = T::lua_constructor(lua_ctx);
        let actions_constructor = Actions::lua_constructor(lua_ctx);
        let action_constructor = Action::lua_constructor(lua_ctx);

        globals.set("RealVec2", realvec2_constructor).unwrap();
        globals.set("UIntVec2", uintvec2_constructor).unwrap();
        globals.set("IntVec2", intvec2_constructor).unwrap();
        globals.set("Model", user_model_constructor).unwrap();
        globals.set("Actions", actions_constructor).unwrap();
        globals.set("Action", action_constructor).unwrap();



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

        conf.model = globals.get::<_,T>("model").expect("failed reading model parameters");
        conf.actions = globals.get::<_,Actions>("actions").expect("failed reading Actions");

    });
    Ok(conf)
}
