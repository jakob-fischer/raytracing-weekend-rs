use crate::rt_core::*;
use std::sync::Arc;

pub struct Sphere {
    center: Point,
    radius: f64,
    material: Arc<Box<dyn Material + Send + Sync>>,
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
                material: self.material.clone(),
            })
        }
    }
}

impl Sphere {
    pub fn new(center: Point, radius: f64, material: Arc<Box<dyn Material + Send + Sync>>) -> Self {
        Sphere {
            center,
            radius,
            material: material.clone(),
        }
    }
}
