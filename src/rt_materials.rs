use crate::rt_core::*;
use rand::prelude::ThreadRng;
use rand::Rng;
use rand::distributions::Uniform;

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

fn reflectance( cosine : f64, ref_idx : f64) -> f64 {
    // Use Schlick's approximation for reflectance.
    let r0 = (1.0-ref_idx) / (1.0+ref_idx);
    let r0 = r0*r0;
    r0 + (1.0-r0)*f64::powf(1.0 - cosine,5.0)
}


pub struct Dielectric {
    ir: f64,
}

impl Dielectric {
    pub fn new(ir: f64) -> Self {
        Self { ir }
    }
}

impl Material for Dielectric {
    fn scatter(
        &self,
        r_in: &Ray,
        hit_record: &HitRecord,
        rng: &mut ThreadRng,
    ) -> Option<ScatterRecord> {
        let unit_direction = r_in.direction.get_normalized();
        let refraction_ratio = if hit_record.front_face {
            1.0 / self.ir
        } else {
            self.ir
        };
        let cos_theta = (-unit_direction.dot(&hit_record.normal)).min(1.0);
        let sin_theta = (1.0 - cos_theta * cos_theta).sqrt();
        let cannot_refract = refraction_ratio * sin_theta > 1.0;

        let dist = Uniform::new(0.0, 1.0);
        let direction = if cannot_refract || reflectance(cos_theta, refraction_ratio) > rng.sample(dist) {
            unit_direction.reflect(&hit_record.normal)
        } else {
            unit_direction.refract(&hit_record.normal, refraction_ratio)
        };

        let scattered_ray = Ray::new(hit_record.point, direction);

        Some(ScatterRecord {
            attenuation: Colour::new(1.0, 1.0, 1.0),
            scattered_ray,
        })
    }
}
