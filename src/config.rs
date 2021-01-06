use rlua::{Lua};
use std::fs;

/// struct which stores all the global configuration for the simulation
#[derive(Debug)]
pub struct Config {
    pub dim: usize,
}

/// function which executes a lua file (located at lua_loc), and returns the configuration
pub fn read_lua(lua_loc: &str) -> Result<Config, std::io::Error>{
    // load the contents of the file
    let lua_file = fs::read_to_string(lua_loc).expect("Could not read lua file");

    // create the lua struct
    let lua = Lua::new();

    // initialise the configuration struct
    let mut conf = Config{
        dim: 0,
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
        conf.dim = globals.get::<_,usize>("dim")
                          .expect("failed setting dimensionality");
    });
    Ok(conf)
}
