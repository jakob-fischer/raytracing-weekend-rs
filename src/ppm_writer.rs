use std::fs::File;
use std::io::Write;

use super::array::Array2d;
use super::math::*;

pub fn write(image: &Array2d<Vec3<f64>>, filename: &str) {
    let mut f = File::create(filename).expect("Unable to create file");
    write!(
        &mut f,
        "P3\n{} {}\n255\n",
        image.get_width(),
        image.get_height()
    )
    .unwrap();

    let mut result = String::new();

    for y in 0..image.get_height() {
        for x in 0..image.get_width() {
            let pixel = image.get(x, y);

            result += &format!(
                "{} {} {}\n",
                (pixel.t[0] * 255.0) as u8,
                (pixel.t[1] * 255.0) as u8,
                (pixel.t[2] * 255.0) as u8
            );
        }
    }

    write!(&mut f, "{}", &result).unwrap();

    f.flush().unwrap();
}
