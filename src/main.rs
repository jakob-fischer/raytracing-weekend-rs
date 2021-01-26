mod array;
mod math;
mod ppm_writer;

use array::Array2d;
use math::*;

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
}

trait Hittable {
    fn hit(&self, ray: &Ray, r: (f64, f64)) -> Option<HitRecord>;
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

            Some(HitRecord {
                point,
                normal: (point - self.center) * (1.0 / self.radius),
                t: root,
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

    fn ray_color(&self) -> Colour {
        if let Some(record) = Sphere::new(&Point::new(0.0, 0.0, -1.0), 0.5).hit(self, (0.0, 1000.0))
        {
            return (record.normal + Colour::new(1.0, 1.0, 1.0)) * 0.5;
        }

        let unit_direction = self.direction.get_normalized();
        let t = 0.5 * (unit_direction.dot(&Direction::new(0.0, 1.0, 0.0)) + 1.0);
        Colour::new(1.0, 1.0, 1.0) * (1.0 - t) + Colour::new(0.5, 0.7, 1.0) * t
    }
}

fn main() {
    // image
    let aspect_ratio = 16.0 / 8.0;
    let image_width = 400 as usize;
    let image_height = (image_width as f64 / aspect_ratio) as usize;

    // camera
    let viewport_height = 2.0;
    let viewport_width = aspect_ratio * viewport_height;
    let focal_length = 1.0;

    let origin = Point::new(0.0, 0.0, 0.0);
    let horizontal = Direction::new(viewport_width, 0.0, 0.0);
    let vertical = Direction::new(0.0, viewport_height, 0.0);
    let lower_left_corner =
        origin - horizontal * 0.5 - vertical * 0.5 - Point::new(0.0, 0.0, focal_length);

    let mut array: Array2d<Colour> =
        Array2d::new(image_width, image_height, &Colour::new(0.0, 0.0, 0.0));
    for x in 0..image_width {
        print!(".");
        for y in 0..image_height {
            let u = x as f64 / (image_width - 1) as f64;
            let v = y as f64 / (image_height - 1) as f64;
            let r = Ray::new(
                origin,
                lower_left_corner + horizontal * u + vertical * v - origin,
            );
            *array.get_mut(x, y) = r.ray_color();
        }
    }
    println!("");

    ppm_writer::write(&array, "test.ppm");
}
