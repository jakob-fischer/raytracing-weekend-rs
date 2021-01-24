use std::ops;
use std::ops::*;

#[derive(Copy, Clone, Debug)]
pub struct Vec3<T: Copy> {
    pub t: [T; 3],
}

impl<T: Copy> Vec3<T> {
    pub fn new(x: T, y: T, z: T) -> Self {
        Self { t: [x, y, z] }
    }
}

impl<T: Add<T, Output = T> + Mul<T, Output = T> + Clone + Copy> Vec3<T> {
    fn squared_length(self: &Vec3<T>) -> T {
        self.t[0] * self.t[0] + self.t[1] * self.t[1] + self.t[2] * self.t[2]
    }
}

pub trait Norm {
    type Length;
    fn length(self: &Self) -> Self::Length;
}

impl Norm for Vec3<f32> {
    type Length = f32;
    fn length(self: &Self) -> Self::Length {
        self.squared_length().sqrt()
    }
}

impl Norm for Vec3<f64> {
    type Length = f64;
    fn length(self: &Self) -> Self::Length {
        self.squared_length().sqrt()
    }
}

impl<T: ops::Neg<Output = T> + Copy> ops::Neg for &Vec3<T> {
    type Output = Vec3<T>;

    fn neg(self) -> Self::Output {
        Vec3::<T> {
            t: [-self.t[0], -self.t[1], -self.t[2]],
        }
    }
}

impl<T: ops::Add<T, Output = T> + Copy> ops::Add<&Vec3<T>> for &Vec3<T> {
    type Output = Vec3<T>;

    fn add(self, rhs: &Vec3<T>) -> Self::Output {
        Vec3::<T> {
            t: [
                self.t[0] + rhs.t[0],
                self.t[1] + rhs.t[1],
                self.t[2] + rhs.t[2],
            ],
        }
    }
}

impl<T: ops::Add<T, Output = T> + Copy> ops::Add<Vec3<T>> for Vec3<T> {
    type Output = Vec3<T>;

    fn add(mut self, rhs: Vec3<T>) -> Self::Output {
        self.t[0] = self.t[0] + rhs.t[0];
        self.t[1] = self.t[1] + rhs.t[1];
        self.t[2] = self.t[2] + rhs.t[2];
        self
    }
}

impl<T: ops::Add<T, Output = T> + Copy> ops::Add<&Vec3<T>> for Vec3<T> {
    type Output = Vec3<T>;

    fn add(mut self, rhs: &Vec3<T>) -> Self::Output {
        self.t[0] = self.t[0] + rhs.t[0];
        self.t[1] = self.t[1] + rhs.t[1];
        self.t[2] = self.t[2] + rhs.t[2];
        self
    }
}

impl<T: ops::Add<T, Output = T> + Copy> ops::Add<Vec3<T>> for &Vec3<T> {
    type Output = Vec3<T>;

    fn add(self, mut rhs: Vec3<T>) -> Self::Output {
        rhs.t[0] = self.t[0] + rhs.t[0];
        rhs.t[1] = self.t[1] + rhs.t[1];
        rhs.t[2] = self.t[2] + rhs.t[2];
        rhs
    }
}

impl<T: ops::AddAssign<T> + Copy> ops::AddAssign<&Vec3<T>> for Vec3<T> {
    fn add_assign(&mut self, rhs: &Self) {
        self.t[0] += rhs.t[0];
        self.t[1] += rhs.t[1];
        self.t[2] += rhs.t[2];
    }
}

impl<T: ops::Sub<T, Output = T> + Copy> ops::Sub<&Vec3<T>> for &Vec3<T> {
    type Output = Vec3<T>;

    fn sub(self, rhs: &Vec3<T>) -> Self::Output {
        Vec3::<T> {
            t: [
                self.t[0] - rhs.t[0],
                self.t[1] - rhs.t[1],
                self.t[2] - rhs.t[2],
            ],
        }
    }
}

impl<T: ops::Sub<T, Output = T> + Copy> ops::Sub<Vec3<T>> for Vec3<T> {
    type Output = Vec3<T>;

    fn sub(mut self, rhs: Vec3<T>) -> Self::Output {
        self.t[0] = self.t[0] - rhs.t[0];
        self.t[1] = self.t[1] - rhs.t[1];
        self.t[2] = self.t[2] - rhs.t[2];
        self
    }
}

impl<T: ops::Sub<T, Output = T> + Copy> ops::Sub<&Vec3<T>> for Vec3<T> {
    type Output = Vec3<T>;

    fn sub(mut self, rhs: &Vec3<T>) -> Self::Output {
        self.t[0] = self.t[0] - rhs.t[0];
        self.t[1] = self.t[1] - rhs.t[1];
        self.t[2] = self.t[2] - rhs.t[2];
        self
    }
}

impl<T: ops::Sub<T, Output = T> + Copy> ops::Sub<Vec3<T>> for &Vec3<T> {
    type Output = Vec3<T>;

    fn sub(self, mut rhs: Vec3<T>) -> Self::Output {
        rhs.t[0] = self.t[0] - rhs.t[0];
        rhs.t[1] = self.t[1] - rhs.t[1];
        rhs.t[2] = self.t[2] - rhs.t[2];
        rhs
    }
}

impl<T: ops::SubAssign<T> + Copy> ops::SubAssign<&Vec3<T>> for Vec3<T> {
    fn sub_assign(&mut self, rhs: &Self) {
        self.t[0] -= rhs.t[0];
        self.t[1] -= rhs.t[1];
        self.t[2] -= rhs.t[2];
    }
}

impl<T: ops::Mul<T, Output = T> + Copy> ops::Mul<&Vec3<T>> for &Vec3<T> {
    type Output = Vec3<T>;

    fn mul(self, rhs: &Vec3<T>) -> Self::Output {
        Vec3::<T> {
            t: [
                self.t[0] * rhs.t[0],
                self.t[1] * rhs.t[1],
                self.t[2] * rhs.t[2],
            ],
        }
    }
}

impl<T: ops::Mul<T, Output = T> + Copy> ops::Mul<Vec3<T>> for Vec3<T> {
    type Output = Vec3<T>;

    fn mul(mut self, rhs: Vec3<T>) -> Self::Output {
        self.t[0] = self.t[0] * rhs.t[0];
        self.t[1] = self.t[1] * rhs.t[1];
        self.t[2] = self.t[2] * rhs.t[2];
        self
    }
}

impl<T: ops::Mul<T, Output = T> + Copy> ops::Mul<Vec3<T>> for &Vec3<T> {
    type Output = Vec3<T>;

    fn mul(self, mut rhs: Vec3<T>) -> Self::Output {
        rhs.t[0] = self.t[0] * rhs.t[0];
        rhs.t[1] = self.t[1] * rhs.t[1];
        rhs.t[2] = self.t[2] * rhs.t[2];
        rhs
    }
}

impl<T: ops::Mul<T, Output = T> + Copy> ops::Mul<&Vec3<T>> for Vec3<T> {
    type Output = Vec3<T>;

    fn mul(mut self, rhs: &Vec3<T>) -> Self::Output {
        self.t[0] = self.t[0] * rhs.t[0];
        self.t[1] = self.t[1] * rhs.t[1];
        self.t[2] = self.t[2] * rhs.t[2];
        self
    }
}

impl<T: ops::MulAssign<T> + Copy> ops::MulAssign<&Vec3<T>> for Vec3<T> {
    fn mul_assign(&mut self, rhs: &Self) {
        self.t[0] *= rhs.t[0];
        self.t[1] *= rhs.t[1];
        self.t[2] *= rhs.t[2];
    }
}

impl<T: ops::Mul<T, Output = T> + Copy> ops::Mul<T> for &Vec3<T> {
    type Output = Vec3<T>;

    fn mul(self, rhs: T) -> Self::Output {
        Vec3::<T> {
            t: [self.t[0] * rhs, self.t[1] * rhs, self.t[2] * rhs],
        }
    }
}

impl<T: ops::Mul<T, Output = T> + Copy> ops::Mul<T> for Vec3<T> {
    type Output = Vec3<T>;

    fn mul(mut self, rhs: T) -> Self::Output {
        self.t[0] = self.t[0] * rhs;
        self.t[1] = self.t[1] * rhs;
        self.t[2] = self.t[2] * rhs;
        self
    }
}

impl<T: ops::MulAssign<T> + Copy> ops::MulAssign<T> for Vec3<T> {
    fn mul_assign(&mut self, rhs: T) {
        self.t[0] *= rhs;
        self.t[1] *= rhs;
        self.t[2] *= rhs;
    }
}


impl<T: ops::Div<T, Output = T> + Copy> ops::Div<&Vec3<T>> for &Vec3<T> {
    type Output = Vec3<T>;

    fn div(self, rhs: &Vec3<T>) -> Self::Output {
        Vec3::<T> {
            t: [
                self.t[0] / rhs.t[0],
                self.t[1] / rhs.t[1],
                self.t[2] / rhs.t[2],
            ],
        }
    }
}

impl<T: ops::DivAssign<T> + Copy> ops::DivAssign<&Vec3<T>> for Vec3<T> {
    fn div_assign(&mut self, rhs: &Self) {
        self.t[0] /= rhs.t[0];
        self.t[1] /= rhs.t[1];
        self.t[2] /= rhs.t[2];
    }
}

impl<T: ops::Div<T, Output = T> + Copy> ops::Div<T> for &Vec3<T> {
    type Output = Vec3<T>;

    fn div(self, rhs: T) -> Self::Output {
        Vec3::<T> {
            t: [self.t[0] / rhs, self.t[1] / rhs, self.t[2] / rhs],
        }
    }
}

impl<T: ops::DivAssign<T> + Copy> ops::DivAssign<T> for Vec3<T> {
    fn div_assign(&mut self, rhs: T) {
        self.t[0] /= rhs;
        self.t[1] /= rhs;
        self.t[2] /= rhs;
    }
}


impl<T> Vec3<T>
where
    T: Copy,
    Vec3<T>: Norm<Length = T> + DivAssign<T>
{
    pub fn get_normalized(&self) -> Self {
        let mut result = *self;
        result /= self.length();
        result
    }
}

impl<T> Vec3<T>
where
    T: Copy,
    Vec3<T>: Norm<Length = T> + DivAssign<T>,
{
    pub fn normalize(&mut self) {
        let len = self.length();
        self.div_assign(len);
    }
}