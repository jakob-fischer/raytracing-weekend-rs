use maglio::*;
use crate::rt_core::*;
use std::sync::Arc;

use maglio::Point3d as Point;
use maglio::ConstrainedRay3d as ConstrainedRay;

pub struct Sphere {
    center: Point,
    radius: f64,
    material: Arc<MaterialBox>,
}

impl Hittable for Sphere {
    fn hit(&self, cray: &ConstrainedRay) -> Option<HitRecord> {
        let oc = cray.ray.origin - self.center;
        let a = cray.ray.direction.squared_length();
        let half_b = oc.dot(&cray.ray.direction);
        let c = oc.squared_length() - self.radius * self.radius;
        let discriminant = half_b * half_b - a * c;

        if discriminant < 0.0 {
            None
        } else {
            let sqrtd = discriminant.sqrt();
            // Find the nearest root that lies in the acceptable range.
            let mut root = (-half_b - sqrtd) / a;
            if root < cray.range.0 || cray.range.1 < root {
                root = (-half_b + sqrtd) / a;
                if root < cray.range.0 || cray.range.1 < root {
                    return None;
                }
            }

            let point = cray.ray.at(root);
            let outward_normal = (point - self.center) * (1.0 / self.radius);
            let front_face = cray.ray.direction.dot(&outward_normal) < 0.0;

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

    fn get_bounding_box(&self) -> Option<BoundingBox3d> {
        let point = Point::new(self.radius, self.radius, self.radius);
        Some(BoundingBox3d {
            u: self.center - point,
            v: self.center + point,
        })
    }
}

impl Sphere {
    pub fn new(center: Point, radius: f64, material: Arc<MaterialBox>) -> Self {
        Sphere {
            center,
            radius,
            material: material.clone(),
        }
    }
}
