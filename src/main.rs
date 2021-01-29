extern crate rand;
use rand::distributions::Uniform;
use rand::{thread_rng, Rng};
use std::rc::Rc;

mod math;
mod ppm_writer;
mod rt_core;
mod rt_hittables;
mod rt_materials;

use ppm_writer::Array2d;
use rt_core::*;
use rt_hittables::*;
use rt_materials::*;

fn main() {
    // rng
    let mut rng = thread_rng();
    let dist = Uniform::new(0.0, 1.0);

    // image
    let aspect_ratio = 16.0 / 8.0;
    let image_width = 400 as usize;
    let image_height = (image_width as f64 / aspect_ratio) as usize;
    let sample_number = 100;
    let max_depth = 50;

    // world
    let material_ground =
        Rc::<Box<dyn Material>>::new(Box::new(Lambertian::new(&Colour::new(0.8, 0.8, 0.0))));
    let material_center =
        Rc::<Box<dyn Material>>::new(Box::new(Lambertian::new(&Colour::new(0.1, 0.2, 0.5))));

    let material_left = Rc::<Box<dyn Material>>::new(Box::new(Dielectric::new(1.5)));
    let material_right =
        Rc::<Box<dyn Material>>::new(Box::new(Metal::new(&Colour::new(0.8, 0.6, 0.2), 0.0)));

    let mut world = HittableList::new();
    world.add(Box::new(Sphere::new(
        Point::new(0.0, 0.0, -1.0),
        0.5,
        material_center,
    )));
    world.add(Box::new(Sphere::new(
        Point::new(0.0, -100.5, -1.0),
        100.0,
        material_ground,
    )));

    world.add(Box::new(Sphere::new(
        Point::new(-1.0, 0.0, -1.0),
        0.5,
        material_left.clone(),
    )));
    world.add(Box::new(Sphere::new(
        Point::new(-1.0, 0.0, -1.0),
        -0.4,
        material_left,
    )));
    world.add(Box::new(Sphere::new(
        Point::new(1.0, 0.0, -1.0),
        0.5,
        material_right,
    )));

    // camera
    let camera = Camera::new(
        Point::new(-2.0, 2.0, 1.0),
        Point::new(0.0, 0.0, -1.0),
        Point::new(0.0, 1.0, 0.0),
        90.0,
        aspect_ratio,
    );

    let mut array: Array2d<Colour> =
        Array2d::new(image_width, image_height, &Colour::new(0.0, 0.0, 0.0));
    for x in 0..image_width {
        println!("{} lines remaining", image_width - x);
        for y in 0..image_height {
            let colour: Colour = (0..100)
                .map(|_| {
                    let u = (x as f64 + rng.sample(dist)) / (image_width - 1) as f64;
                    let v = (y as f64 + rng.sample(dist)) / (image_height - 1) as f64;
                    camera.get_ray(u, v).ray_color(&world, &mut rng, max_depth)
                })
                .fold(Colour::new(0.0, 0.0, 0.0), |old, n| old + n);

            *array.get_mut(x, y) = colour * (1.0 / sample_number as f64);
        }
    }
    println!("");

    ppm_writer::write(&array, "test.ppm");
}
