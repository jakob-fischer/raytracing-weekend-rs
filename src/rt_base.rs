use crate::math::*;

pub type Colour = Vec3d;
pub type Point = Vec3d;
pub type Direction = Vec3d;

#[derive(Clone)]
pub struct Ray {
    pub origin: Point,
    pub direction: Direction,
}

pub struct ConstrainedRay{
    pub ray : Ray,
    pub range : (f64, f64),
}

#[derive(Clone)]
pub struct BoundingBox3d {
    pub u: Vec3d,
    pub v: Vec3d,
}

impl BoundingBox3d {
    pub fn intersects_in_dimentions(&self, other: &Self, dimension: usize) -> bool {
        self.v.t[dimension] >= other.u.t[dimension] && self.u.t[dimension] <= other.v.t[dimension]
    }

    fn intersects(&self, ray : &ConstrainedRay) -> bool{

        false
    }
}