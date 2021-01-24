use std::fs::File;
use std::io::Write;

mod math;
mod array;

use math::*;
use array::Array2d;


type Colour = Vec3<f32>;
fn write_to_ppm(image: &Array2d<Colour>, filename: &str) {
    let mut f = File::create(filename).expect("Unable to create file");
    write!(
        &mut f,
        "P3\n{} {}\n255\n",
        image.get_width(),
        image.get_height()
    );

    for y in 0..image.get_height() {
        for x in 0..image.get_width() {
            let pixel = image.get(x, y);

            write!(
                &mut f,
                "{} {} {}\n",
                (pixel.t[0] * 255.0) as u8,
                (pixel.t[1] * 255.0) as u8,
                (pixel.t[2] * 255.0) as u8
            ).unwrap();
        }
    }

    f.flush().unwrap();
}

fn main() {
    let nx = 200;
    let ny = 100;
    let mut array: Array2d<Colour> = Array2d::new_init(nx, ny, &Colour::new(0.0, 0.0, 0.0));

    for x in 0..nx {
        for y in 0..ny {
            *array.get_mut(x, y) =
                Colour::new((x as f32) / (nx as f32), (y as f32) / (ny as f32), 0.2);
        }
    }

    println!("Hello, world!");

    write_to_ppm(&array, "test.ppm");
}
