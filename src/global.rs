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
pub struct GlobalVal
{
    pub value: Real,
    get_value: fn(&CartesianBlock)->Real,
}
impl GlobalVal{
    pub fn calc_value(&mut self, dfs: &CartesianBlock)
    {
        self.value = (self.get_value)(dfs);
    }
}


