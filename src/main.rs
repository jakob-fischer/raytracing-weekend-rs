extern crate rand;
use rand::distributions::Uniform;
use rand::prelude::ThreadRng;
use rand::{thread_rng, Rng};
use rayon::prelude::*;
use std::sync::Arc;

mod math;
mod ppm_writer;
mod rt_core;
mod rt_hittables;
mod rt_materials;

use math::*;
use ppm_writer::Array2d;
use rt_core::*;
use rt_hittables::*;
use rt_materials::*;

fn random_double(min: f64, max: f64, rng: &mut ThreadRng) -> f64 {
    let dist = Uniform::new(min, max);
    rng.sample(dist)
}

fn sample(rng: &mut ThreadRng) -> f64 {
    let dist = Uniform::new(0.0, 1.0);
    rng.sample(dist)
}

fn random_vec(min: f64, max: f64, rng: &mut ThreadRng) -> Vec3<f64> {
    let dist = Uniform::new(min, max);
    Vec3::<f64>::new(rng.sample(dist), rng.sample(dist), rng.sample(dist))
}

fn random_scene(rng: &mut ThreadRng) -> Arc<HittableBox> {
    let material1 = Arc::new(Box::new(Dielectric::new(1.5)) as MaterialBox);
    let material2 = Arc::new(Box::new(Lambertian::new(&Colour::new(0.4, 0.2, 0.1))) as MaterialBox);
    let material3 = Arc::new(Box::new(Metal::new(&Colour::new(0.7, 0.6, 0.5), 0.0)) as MaterialBox);
    let ground_material =
        Arc::new(Box::new(Lambertian::new(&Colour::new(0.5, 0.5, 0.5))) as MaterialBox);

    let mut world = HittableList::new();

    world.add(Box::new(Sphere::new(
        Point::new(0.0, -1000.0, 0.0),
        1000.0,
        ground_material,
    )));

    for a in -11..=11 {
        for b in -11..11 {
            let choose_mat = sample(rng);
            let center = Point::new(
                a as f64 + 0.9 * sample(rng),
                0.2,
                b as f64 + 0.9 * sample(rng),
            );

            if (center - Point::new(4.0, 0.2, 0.0)).length() > 0.9 {
                let material = if choose_mat < 0.8 {
                    // diffuse
                    let albedo = random_vec(0.0, 1.0, rng) * random_vec(0.0, 1.0, rng);
                    Arc::new(Box::new(Lambertian::new(&albedo)) as MaterialBox)
                } else if choose_mat < 0.95 {
                    // metal
                    let albedo = random_vec(0.5, 1.0, rng);
                    let fuzz = random_double(0.0, 0.5, rng);
                    Arc::new(Box::new(Metal::new(&albedo, fuzz)) as MaterialBox)
                } else {
                    // glass
                    Arc::new(Box::new(Dielectric::new(1.5)) as MaterialBox)
                };

                world.add(Box::new(Sphere::new(center, 0.2, material)));
            }
        }
    }

    world.add(Box::new(Sphere::new(Point::new(0.0, 1.0, 0.0), 1.0, material1)) as HittableBox);

    world.add(Box::new(Sphere::new(
        Point::new(-4.0, 1.0, 0.0),
        1.0,
        material2,
    )));

    world.add(Box::new(Sphere::new(
        Point::new(4.0, 1.0, 0.0),
        1.0,
        material3,
    )));
    Arc::new(Box::new(world) as HittableBox)
}

fn main() {
    // rng
    let mut rng = thread_rng();
    let dist = Uniform::new(0.0, 1.0);

    // image
    let aspect_ratio = 3.0 / 2.0;
    let image_width = 1200;
    let image_height = (image_width as f64 / aspect_ratio) as usize;
    let sample_number = 5;
    let max_depth = 50;

    // world
    let world = random_scene(&mut rng);

    // camera
    let lookfrom = Point::new(13.0, 2.0, 3.0);
    let lookat = Point::new(0.0, 0.0, 0.0);
    let vup = Point::new(0.0, 1.0, 0.0);
    let dist_to_focus = 10.0;
    let aperture = 0.1;

    let camera = Camera::new(
        lookfrom,
        lookat,
        vup,
        20.0,
        aspect_ratio,
        aperture,
        dist_to_focus,
    );

    let mut array: Array2d<Colour> =
        Array2d::new(image_width, image_height, &Colour::new(0.0, 0.0, 0.0));

    let pixel_indexes: Vec<_> = (0..image_width)
        .rev()
        .map(|x| -> Vec<_> { (0..image_height).map(|y| (x, y)).collect() })
        .flatten()
        .collect();

    let result : Vec<_> = pixel_indexes.par_iter().map(|(x, y)| {
        let mut rng = thread_rng();

        let colour: Colour = (0..sample_number)
            .map(|_| {
                let u = (*x as f64 + rng.sample(dist)) / (image_width - 1) as f64;
                let v = (*y as f64 + rng.sample(dist)) / (image_height - 1) as f64;
                camera
                    .get_ray(u, v, &mut rng)
                    .ray_color(world.clone(), &mut rng, max_depth)
            })
            .fold(Colour::new(0.0, 0.0, 0.0), |old, n| old + n);

        (*x, *y, colour * (1.0 / sample_number as f64))
    }).collect();

    for (x, y, colour) in result {
        *array.get_mut(x, y) = colour;
    }

    ppm_writer::write(&array, "test.ppm");
}
