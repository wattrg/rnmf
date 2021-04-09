use rlua::{Lua, UserData, Context};
use std::fs;
use crate::*;
use colored::*;
use std::process::Command;
//use crate::solver::Actions;
use crate::io::OutputCallBackHashMap;

pub trait UserConfig: UserData + Clone {
    fn new()->Self;
    fn lua_constructor(lua_ctx: Context)->rlua::Function;
}


pub struct Config<T: UserConfig>{
    pub geom: GeomConfig,
    pub model: T,
    //pub actions: Actions<'a, T>,
    pub residual_iters: usize,
}

impl <T: UserConfig> Config<T>{
    pub fn new(user_config: T)->Config<T>{
        Config{
            //geom: GeomConfig::new(),
            geom: Default::default(),
            model: user_config,
            //actions: Actions::new(),
            residual_iters: 1,
        }
    }
}

#[derive(Default)]
pub struct GeomConfig{
    pub dim: usize,
    pub length: RealVec2,
    pub n_cells: UIntVec2,
}
// impl GeomConfig{
//     pub fn new() -> GeomConfig{
//         GeomConfig{
//             dim: 0,
//             length: RealVec2([0.0, 0.0]),
//             n_cells: UIntVec2([0, 0]),
//         }
//     }
// }

impl UserData for RealVec2 {}
impl UserData for UIntVec2 {}
impl UserData for IntVec2 {}
impl UserData for RealVec3 {}
impl UserData for UIntVec3 {}
impl UserData for IntVec3 {}

pub fn init<T: 'static>(args: Vec<String>, user_model: T, out_cb:fn()->OutputCallBackHashMap)
    ->Result<(Config<T>, OutputCallBackHashMap), std::io::Error>
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


    if args.len() <= 1 {
        println!("{} Lua configuration script not given.", "Error:".red());
        panic!();
    }
    

    Ok((read_lua(args[1].clone(), user_model)?, out_cb()))


}

/// function which executes a lua file (located at lua_loc), and returns the configuration
pub fn read_lua<T: 'static>(lua_loc: String, user_model: T) -> Result<Config<T>, std::io::Error>
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
        //let actions_constructor = Actions::lua_constructor(lua_ctx);
        //let action_constructor = Action::lua_constructor(lua_ctx);

        globals.set("RealVec2", realvec2_constructor).unwrap();
        globals.set("UIntVec2", uintvec2_constructor).unwrap();
        globals.set("IntVec2", intvec2_constructor).unwrap();
        globals.set("Model", user_model_constructor).unwrap();
        //globals.set("Actions", actions_constructor).unwrap();
        //globals.set("Action", action_constructor).unwrap();



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
        //conf.actions = globals.get::<_,Actions>("actions").expect("failed reading Actions");

    });
    Ok(conf)
}
