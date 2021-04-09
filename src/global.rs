use crate::Real;
use crate::mesh::cartesian2d::CartesianBlock;

#[derive(Clone)]
pub enum GlobalCompare {
    Greater(Real),
    Less(Real),
    GreaterEq(Real),
    LessEq(Real),
    Eq(Real),
}

#[derive(Clone)]
pub struct GlobalVal<T>
{
    pub value: Real,
    get_value: fn(&CartesianBlock<T>)->Real,
}
impl <T> GlobalVal<T>{
    pub fn calc_value(&mut self, dfs: &CartesianBlock<T>)
    {
        self.value = (self.get_value)(dfs);
    }
}


