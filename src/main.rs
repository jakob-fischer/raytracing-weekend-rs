extern crate rand;
use rand::distributions::Uniform;
use rand::{thread_rng, Rng, RngCore};

mod array;
mod math;
mod ppm_writer;

use array::Array2d;
use math::*;

fn degrees_to_radians(degrees: f64) -> f64 {
    return degrees * std::f64::consts::PI / 180.0;
}

type Colour = Vec3<f64>;
type Point = Vec3<f64>;
type Direction = Vec3<f64>;

struct Ray {
    origin: Point,
    direction: Direction,
}

struct HitRecord {
    point: Point,
    normal: Direction,
    t: f64,
    front_face: bool,
}

trait Hittable {
    fn hit(&self, ray: &Ray, r: (f64, f64)) -> Option<HitRecord>;
}

struct HittableList {
    hittable: Vec<Box<dyn Hittable>>,
}

impl HittableList {
    fn new() -> HittableList {
        HittableList {
            hittable: Vec::new(),
        }
    }

    fn clear(&mut self) {
        self.hittable.clear();
    }

    fn add(&mut self, hittable: Box<dyn Hittable>) {
        self.hittable.push(hittable);
    }
}

impl Hittable for HittableList {
    fn hit(&self, ray: &Ray, r: (f64, f64)) -> Option<HitRecord> {
        self.hittable
            .iter()
            .map(|x| x.hit(ray, (0.0, f64::INFINITY)))
            .filter(|x| x.is_some())
            .map(|x| x.unwrap())
            .min_by(|x, y| x.t.partial_cmp(&y.t).unwrap())
    }
}

struct Sphere {
    center: Point,
    radius: f64,
}

impl Hittable for Sphere {
    fn hit(&self, ray: &Ray, (t_min, t_max): (f64, f64)) -> Option<HitRecord> {
        let oc = ray.origin - self.center;
        let a = ray.direction.squared_length();
        let half_b = oc.dot(&ray.direction);
        let c = oc.squared_length() - self.radius * self.radius;
        let discriminant = half_b * half_b - a * c;

        if discriminant < 0.0 {
            None
        } else {
            let sqrtd = discriminant.sqrt();
            // Find the nearest root that lies in the acceptable range.
            let mut root = (-half_b - sqrtd) / a;
            if root < t_min || t_max < root {
                root = (-half_b + sqrtd) / a;
                if root < t_min || t_max < root {
                    return None;
                }
            }

            let point = ray.at(root);
            let outward_normal = (point - self.center) * (1.0 / self.radius);
            let front_face = ray.direction.dot(&outward_normal) < 0.0;

            Some(HitRecord {
                point,
                normal: if front_face {
                    outward_normal
                } else {
                    -outward_normal
                },
                t: root,
                front_face,
            })
        }
    }
}

impl Sphere {
    fn new(center: &Point, radius: f64) -> Self {
        Sphere {
            center: *center,
            radius,
        }
    }
}

impl Ray {
    fn new(origin: Point, direction: Direction) -> Self {
        Self { origin, direction }
    }

    fn at(&self, t: f64) -> Point {
        &(&self.direction * t) + &self.origin
    }

    fn ray_color(&self, world: &dyn Hittable) -> Colour {
        if let Some(record) = world.hit(self, (0.0, f64::INFINITY)) {
            return (record.normal + Colour::new(1.0, 1.0, 1.0)) * 0.5;
        }

        let unit_direction = self.direction.get_normalized();
        let t = 0.5 * (unit_direction.dot(&Direction::new(0.0, 1.0, 0.0)) + 1.0);
        Colour::new(1.0, 1.0, 1.0) * (1.0 - t) + Colour::new(0.5, 0.7, 1.0) * t
    }
}

struct Camera {
    viewport_height: f64,
    viewport_width: f64,
    focal_length: f64,
    origin: Point,
    horizontal: Direction,
    vertical: Direction,
    lower_left_corner: Vec3<f64>,
}

impl Camera {
    fn new(aspect_ratio: f64) -> Self {
        let viewport_height = 2.0;
        let viewport_width = aspect_ratio * viewport_height;
        let origin = Point::new(0.0, 0.0, 0.0);
        let horizontal = Direction::new(viewport_width, 0.0, 0.0);
        let vertical = Direction::new(0.0, viewport_height, 0.0);
        let focal_length = 1.0;
        let lower_left_corner =
            origin - horizontal * 0.5 - vertical * 0.5 - Point::new(0.0, 0.0, focal_length);

        Camera {
            viewport_height,
            viewport_width,
            focal_length,
            origin,
            horizontal,
            vertical,
            lower_left_corner,
        }
    }

    fn get_ray(&self, u: f64, v: f64) -> Ray {
        Ray::new(
            self.origin,
            self.lower_left_corner + self.horizontal * u + self.vertical * v - self.origin,
        )
    }
}

fn main() {
    // rng
    let mut rng = thread_rng();
    let dist = Uniform::new(0.0, 1.0);

    // image
    let aspect_ratio = 16.0 / 8.0;
    let image_width = 400 as usize;
    let image_height = (image_width as f64 / aspect_ratio) as usize;
    let sample_number = 100;

    // world
    let mut world = HittableList::new();
    world.add(Box::new(Sphere::new(&Point::new(0.0, 0.0, -1.0), 0.5)));
    world.add(Box::new(Sphere::new(&Point::new(0.0, -100.5, -1.0), 100.0)));

    // camera
    let camera = Camera::new(aspect_ratio);

    let mut array: Array2d<Colour> =
        Array2d::new(image_width, image_height, &Colour::new(0.0, 0.0, 0.0));
    for x in 0..image_width {
        print!(".");
        for y in 0..image_height {
            let colour: Colour = (0..100)
                .map(|_| {
                    let u = (x as f64 + rng.sample(dist)) / (image_width - 1) as f64;
                    let v = (y as f64 + rng.sample(dist)) / (image_height - 1) as f64;
                    camera.get_ray(u, v).ray_color(&world)
                })
                .fold(Colour::new(0.0, 0.0, 0.0), |old, n| old + n);

            *array.get_mut(x, y) = colour * (1.0 / sample_number as f64);
        }
    }
    println!("");

    ppm_writer::write(&array, "test.ppm");
}
