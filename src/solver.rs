use crate::global::{GlobalCompare, GlobalVal};
use crate::config::{Config, UserConfig};
use crate::mesh::cartesian2d::CartesianBlock;


/// A stop condition tells the solver when to stop iterating
/// There must be a maximum number of iterations, `iters`, and optionally can have 
/// a stop condition defined based on a global quantity. 
pub struct Stop<'a>
{
    pub max_iters: usize,
    pub global_stop: Option<&'a mut [GlobalStop]>,
}
impl <'a> Stop<'a> {
    fn new() -> Stop<'a> {
        Stop{
            max_iters: 0,
            global_stop: None,
        }
    }
}

/// `GlobalStop` contains the mechanism to check if the stop condition has been met
#[derive(Clone)]
pub struct GlobalStop
{
    pub global: GlobalVal,
    pub condition: GlobalCompare,
    pub check_every: usize,
}
impl GlobalStop {
    fn check_stop_condition(&self, compare_type: &GlobalCompare) -> bool {
        match compare_type {
            GlobalCompare::Greater(threshold) => { self.global.value > *threshold }
            GlobalCompare::Less(threshold) => { self.global.value < *threshold }
            GlobalCompare::GreaterEq(threshold) => { self.global.value >= *threshold }
            GlobalCompare::LessEq(threshold) => { self.global.value <= *threshold }
            GlobalCompare::Eq(threshold) => { self.global.value == *threshold }
        }
    }
}

/// `Actions` stores a list of `Action` structs. Each `Action` is executed sequentially.
/// An example of using multiple actions the first action might progress a field to steady state, 
/// and then the second action will propagate both fields forwards at the same time. This allows 
/// for run-time configuration of what actually gets executed.
pub struct Actions<'a, T:UserConfig> (Vec<Action<'a, T>>);
impl <'a, T> Actions<'a, T> 
    where T: UserConfig
{
    pub fn new()->Actions<'a,T>{
        Actions(Vec::new())
    }
}
impl <'a, T> core::ops::Deref for Actions<'a, T>
    where T: UserConfig
{
    type Target = Vec<Action<'a, T>>;

    fn deref (self: &'_ Self) -> &'_ Self::Target {
        &self.0
    }
}

/// An `Action` contains the information to complete one of the steps in the `Actions`. 
/// Each `Action` has three elements (or at least they will)
/// - `pre_action` is run once before the main action starts. A possible use might be to modify 
/// a setting in the configuration. For example, a computational fluid dynamicist may wish 
/// to run one action with viscosity turn off until the flow reaches steady state, 
/// then turn viscosity and continue the simulation. 
/// This could be done by using two `Action`s, and having `pre_action` on the second `Action`
///  turn viscosity on.
/// - `action` is the primary iteration `Stage`. It contains a sequence of `Step`s, 
/// which will be executed one after the other. 
/// Each `Stage` of each `Step` contain another `Stage`. So if you need an iterative process 
/// within an iterative process within an iterative process, its turtles all the way down.
/// Each `Stage` contains a stopping condition, which can depend on global variables meeting 
/// some condition.
/// - `post_action` is something which is executed immediately after the main iteration 
/// in the `action`
pub struct Action<'a, T: UserConfig> {
    // pre_action: ,

    /// Stores how to complete one iteration, and when to stop iterating
    action: Stage<'a, T>,

    // post_action: ,
}
impl <'a, T> Action<'a, T>
    where T: UserConfig
{
    pub fn new() -> Action<'a, T> {
        Action{
            action: Stage::new(),
        }
    }
}

/// A stage contains the sequence of `Step`s to complete (called `Iteration`), and the stop condition
pub struct Stage<'a, T: UserConfig> 
{
    /// A sequence of `Steps`
    iteration: Iteration<'a, T>,

    /// The stopping condition
    stop: Stop<'a>,
}
impl <'a, T> Stage<'a, T> 
    where T: UserConfig
{
    pub fn new() -> Stage<'a, T> {
        Stage{
            iteration: Iteration::new(),
            stop: Stop::new(),
        }
    }

    fn execute(&mut self, config: &Config<T>, data_buf: &mut CartesianBlock) {
        self.iteration.execute(config, data_buf);
    }
}

/// Stores the sequence of `Step`s to complete one iteration (which may involve iterative processes)
pub struct Iteration<'a, T: UserConfig>(Vec<Step<'a, T>>);
impl <'a, T> core::ops::Deref for Iteration<'a, T>
    where T: UserConfig
{
    type Target = Vec<Step<'a, T>>;

    fn deref (self: &'_ Self) -> &'_ Self::Target {
        &self.0
    }
}
impl <'a, T> core::ops::DerefMut for Iteration<'a,T>
    where T: UserConfig
{
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}
impl <'a, T> Iteration<'a, T> 
    where T: UserConfig
{
    fn new() -> Iteration<'a, T> {
        Iteration(Vec::new())
    }

    fn execute(&mut self, config: &Config<T>, data_buf: &mut CartesianBlock){
        let mut iter = 0;
        for step in self.iter_mut(){
            iter += 1;
            step.execute(config, data_buf);

            if iter % config.residual_iters == 0 {
                // compute all the globals and residuals
            }
        }
    }
}

/// A `Step` can either be a simple one step calculation (`Callback`) or it can be another `Stage`
pub enum Step<'a, T: UserConfig>
{
    /// A calculation which doesn't require iteration
    Callback( fn(&Config<T>, &mut CartesianBlock) ),

    /// A calculation which requires iteration
    Stage(Stage<'a, T>),
}

impl <'a,T> Step<'a,T>
    where T: UserConfig
{
    fn execute(& mut self, config: &Config<T>, data_buf: &mut CartesianBlock){
        match self{
            // the step consists of a single step... nice and easy!
            Step::Callback(call_back) => {call_back(config, data_buf);}

            // the step requires iteration... grrrrr
            Step::Stage(stage) => {
                for iter_num in 0..stage.stop.max_iters {
                    // execute the stage
                    stage.execute(config, data_buf);

                    // check if we should stop
                    match &mut stage.stop.global_stop{
                        Some(stops) => {
                            for stop in stops.iter_mut(){
                                if stop.check_every % iter_num == 0 {
                                    stop.global.calc_value(data_buf);
                                    if stop.check_stop_condition(&stop.condition){
                                        break;
                                    }
                                }
                            }  
                        }
                        None => {}
                    }
                }
            }
        }
    }
}
