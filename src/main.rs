mod math;
mod array;
mod ppm_writer;

use math::*;
use array::Array2d;


type Colour = Vec3<f64>;
type Point = Vec3<f64>;

fn main() {
    let nx = 200;
    let ny = 100;
    let mut array: Array2d<Colour> = Array2d::new_init(nx, ny, &Colour::new(0.0, 0.0, 0.0));

    for x in 0..nx {
        for y in 0..ny {
            *array.get_mut(x, y) =
                Colour::new((x as f64) / (nx as f64), (y as f64) / (ny as f64), 0.2);
        }
    }

    println!("Hello, world!");

    ppm_writer::write(&array, "test.ppm");
}
