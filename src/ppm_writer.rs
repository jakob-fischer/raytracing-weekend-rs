use std::fs::File;
use std::io::Write;

use super::math::*;

pub struct Array2d<T> {
    width: usize,
    height: usize,
    vec: Vec<T>,
}

impl<T: Copy> Array2d<T> {
    pub fn new(width: usize, height: usize, initial: &T) -> Self {
        Array2d::<T> {
            width,
            height,
            vec: (0..(width * height)).map(|_| *initial).collect(),
        }
    }
}

impl<T> Array2d<T> {
    pub fn get_width(&self) -> usize {
        self.width
    }

    pub fn get_height(&self) -> usize {
        self.height
    }

    pub fn get(&self, x: usize, y: usize) -> &T {
        self.vec.get(x * self.height + y).unwrap()
    }

    pub fn get_mut(&mut self, x: usize, y: usize) -> &mut T {
        self.vec.get_mut(x * self.height + y).unwrap()
    }
}

fn clamp(input: f64, min: f64, max: f64) -> f64 {
    if input < min {
        min
    } else if input > max {
        max
    } else {
        input
    }
}

pub fn write(image: &Array2d<Vec3<f64>>, filename: &str) {
    let mut f = File::create(filename).expect("Unable to create file");
    write!(
        &mut f,
        "P3\n{} {}\n255\n",
        image.get_width(),
        image.get_height()
    )
    .unwrap();

    let get_string_for_pixel = |(x, y)| {
        let pixel = image.get(x, y);
        let clamp = |x: f64| clamp(x.sqrt(), 0.0, 0.999);

        format!(
            "{} {} {}\n",
            (clamp(pixel.t[0]) * 255.0) as u8,
            (clamp(pixel.t[1]) * 255.0) as u8,
            (clamp(pixel.t[2]) * 255.0) as u8
        )
    };

    let result = (0..image.get_height())
        .rev()
        .map(|y| -> Vec<_> { (0..image.get_width()).map(|x| (x, y)).collect() })
        .flatten()
        .map(get_string_for_pixel)
        .fold(String::new(), |acc, next| acc + &next);

    write!(&mut f, "{}", &result).unwrap();
    f.flush().unwrap();
}
