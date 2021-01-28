use crate::rt_core::*;
use rand::prelude::ThreadRng;

pub struct Lambertian {
    albedo: Colour,
}

impl Lambertian {
    pub fn new(albedo: &Colour) -> Self {
        Self { albedo: *albedo }
    }
}

impl Material for Lambertian {
    fn scatter(
        &self,
        _r_in: &Ray,
        hit_record: &HitRecord,
        rng: &mut ThreadRng,
    ) -> Option<ScatterRecord> {
        let scatter_direction = hit_record.normal + random_in_hemisphere(&hit_record.normal, rng);
        let scatter_direction = if scatter_direction.is_almost_zero() {
            hit_record.normal
        } else {
            scatter_direction
        };

        let scattered_ray = Ray::new(hit_record.point, scatter_direction);
        Some(ScatterRecord {
            attenuation: self.albedo,
            scattered_ray,
        })
    }
}

pub struct Metal {
    albedo: Colour,
    fuzz: f64,
}

impl Metal {
    pub fn new(albedo: &Colour, fuzz: f64) -> Self {
        Self {
            albedo: *albedo,
            fuzz,
        }
    }
}

impl Material for Metal {
    fn scatter(
        &self,
        r_in: &Ray,
        hit_record: &HitRecord,
        rng: &mut ThreadRng,
    ) -> Option<ScatterRecord> {
        let reflected = r_in.direction.reflect(&hit_record.normal);
        let scattered_ray = Ray::new(
            hit_record.point,
            reflected + random_unit_sphere(rng) * self.fuzz,
        );
        Some(ScatterRecord {
            attenuation: self.albedo,
            scattered_ray,
        })
    }
}
