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

    let get_string_for_pixel = |(x, y)| {
        let pixel = image.get(x, y);

        format!(
            "{} {} {}\n",
            (pixel.t[0] * 255.0) as u8,
            (pixel.t[1] * 255.0) as u8,
            (pixel.t[2] * 255.0) as u8
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
