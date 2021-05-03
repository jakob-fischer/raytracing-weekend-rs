use crate::kd_tree::*;
use crate::math::*;
use rand::distributions::Uniform;
use rand::prelude::ThreadRng;
use rand::Rng;
use std::sync::Arc;

pub use crate::rt_base::*;

fn degrees_to_radians(degrees: f64) -> f64 {
    return degrees * std::f64::consts::PI / 180.0;
}

impl Vec3d {
    pub fn is_almost_zero(&self) -> bool {
        self.squared_length() < 0e-14
    }

    pub fn reflect(&self, n: &Vec3d) -> Vec3d {
        self - *n * self.dot(n) * 2.0
    }

    pub fn refract(&self, n: &Vec3d, etai_over_etat: f64) -> Vec3d {
        let cos_theta = -n.dot(self).min(1.0);
        let r_out_perp = (self + n * cos_theta) * etai_over_etat;
        let r_out_parallel = n * (-(1.0 - r_out_perp.squared_length()).abs().sqrt());
        r_out_perp + r_out_parallel
    }
}

pub fn random_unit_sphere(rng: &mut ThreadRng) -> Direction {
    let dist = Uniform::new(-1.0, 1.0);
    let mut sample = || rng.sample(dist);

    loop {
        let candidate = Direction::new(sample(), sample(), sample());

        if candidate.squared_length() <= 1.0 {
            return candidate;
        }
    }
}

pub fn random_in_unit_disk(rng: &mut ThreadRng) -> Direction {
    let dist = Uniform::new(-1.0, 1.0);
    let mut sample = || rng.sample(dist);

    loop {
        let candidate = Direction::new(sample(), sample(), 0.0);

        if candidate.squared_length() <= 1.0 {
            return candidate;
        }
    }
}

pub fn random_in_hemisphere(normal: &Direction, rng: &mut ThreadRng) -> Direction {
    let random_in_unit_sphere = random_unit_sphere(rng);

    if normal.dot(&random_in_unit_sphere) > 0.0 {
        random_in_unit_sphere
    } else {
        -random_in_unit_sphere
    }
}

pub struct ScatterRecord {
    pub attenuation: Colour,
    pub scattered_ray: Ray,
}

pub struct HitRecord {
    pub point: Point,
    pub normal: Direction,
    pub t: f64,
    pub front_face: bool,
    pub material: Arc<MaterialBox>,
}

pub trait Hittable {
    fn hit(&self, ray: &ConstrainedRay) -> Option<HitRecord>;
}

pub type HittableBox = Box<dyn Hittable + Send + Sync>;

struct HittableBoxInTree(HittableBox);
type RtKdTree = KdTree<BoundingBox3d, HittableBoxInTree>;

impl KdTreeContent<BoundingBox3d> for HittableBoxInTree {
    fn get_bounding_box(&self) -> BoundingBox3d {
        BoundingBox3d::default()
    }
}

pub struct HittableList {
    hittable: Vec<HittableBox>,
    hittable2: RtKdTree,
}

impl HittableList {
    pub fn new() -> HittableList {
        HittableList {
            hittable: Vec::new(),
            hittable2: RtKdTree::new(),
        }
    }

    pub fn add(&mut self, hittable: HittableBox) {
        self.hittable.push(hittable);
    }
}

impl Hittable for HittableList {
    fn hit(&self, cray: &ConstrainedRay) -> Option<HitRecord> {
        self.hittable
            .iter()
            .map(|x| x.hit(cray))
            .filter(|x| x.is_some())
            .map(|x| x.unwrap())
            .min_by(|x, y| x.t.partial_cmp(&y.t).unwrap())
    }
}

pub trait Material {
    fn scatter(&self, r_in: &Ray, hit: &HitRecord, rng: &mut ThreadRng) -> Option<ScatterRecord>;
}

pub type MaterialBox = Box<dyn Material + Send + Sync>;

impl Ray {
    pub fn new(origin: Point, direction: Direction) -> Self {
        Self { origin, direction }
    }

    pub fn at(&self, t: f64) -> Point {
        &(&self.direction * t) + &self.origin
    }

    pub fn ray_color(&self, world: Arc<HittableBox>, rng: &mut ThreadRng, depth: i32) -> Colour {
        // If we've exceeded the ray bounce limit, no more light is gathered.
        if depth <= 0 {
            Colour::new(0.0, 0.0, 0.0)
        } else if let Some(record) = world.hit(&ConstrainedRay {
            ray: self.clone(),
            range: (0.001, f64::INFINITY),
        }) {
            if let Some(scatter_record) = record.material.scatter(self, &record, rng) {
                scatter_record.attenuation
                    * scatter_record
                        .scattered_ray
                        .ray_color(world, rng, depth - 1)
            } else {
                Colour::new(0.0, 0.0, 0.0)
            }
        } else {
            let unit_direction = self.direction.get_normalized();
            let t = 0.5 * (unit_direction.dot(&Direction::new(0.0, 1.0, 0.0)) + 1.0);
            Colour::new(1.0, 1.0, 1.0) * (1.0 - t) + Colour::new(0.5, 0.7, 1.0) * t
        }
    }
}

pub struct Camera {
    origin: Point,
    horizontal: Direction,
    vertical: Direction,
    lower_left_corner: Vec3d,
    lens_radius: f64,
    u: Vec3d,
    v: Vec3d,
}

impl Camera {
    pub fn new(
        lookfrom: Point,
        lookat: Point,
        vup: Direction,
        vfov: f64,
        aspect_ratio: f64,
        aperture: f64,
        focus_dist: f64,
    ) -> Self {
        let theta = degrees_to_radians(vfov);
        let h = (theta / 2.0).tan();
        let viewport_height = 2.0 * h;
        let viewport_width = aspect_ratio * viewport_height;

        let w = (lookfrom - lookat).get_normalized();
        let u = vup.cross(&w).get_normalized();
        let v = w.cross(&u);

        let origin = lookfrom;
        let horizontal = u * (focus_dist * viewport_width);
        let vertical = v * (focus_dist * viewport_height);
        let lower_left_corner = origin - horizontal * 0.5 - vertical * 0.5 - w * focus_dist;

        let lens_radius = aperture / 2.0;

        Camera {
            origin,
            horizontal,
            vertical,
            lower_left_corner,
            lens_radius,
            u,
            v,
        }
    }

    pub fn get_ray(&self, s: f64, t: f64, rng: &mut ThreadRng) -> Ray {
        let rd = random_in_unit_disk(rng) * self.lens_radius;
        let offset = self.u * rd.t[0] + self.v * rd.t[1];

        Ray::new(
            self.origin + offset,
            self.lower_left_corner + self.horizontal * s + self.vertical * t - self.origin - offset,
        )
    }
}
