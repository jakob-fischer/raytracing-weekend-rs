use crate::math::*;
use crate::rt_base::*;
use std::collections::HashSet;
use std::hash::{Hash, Hasher};
use std::sync::Arc;

pub trait BoundingBoxTrait: Sized + Clone + Default {
    fn partition(&self) -> Option<(Self, Self)>;
    fn extend(&self) -> Option<(Self, Self)>; // First is other, Second new parent

    fn is_contained(&self, other: &Self) -> bool;
    fn intersects(&self, other: &Self) -> bool;

    fn is_sub_scale(&self, other: &Self) -> bool;
}

pub trait HittableBoundingBoxTrait: BoundingBoxTrait {
    fn hit(&self, ray: &Ray, r: (f64, f64)) -> bool;
}

impl Default for BoundingBox3d {
    fn default() -> Self {
        BoundingBox3d {
            u: Vec3d::new(0.0, 0.0, 0.0),
            v: Vec3d::new(1.0, 1.0, 1.0),
        }
    }
}

impl BoundingBox3d {
    fn get_most_narrow_dimension(&self) -> usize {
        let dx = self.u.t[0] - self.v.t[0];
        let dy = self.u.t[1] - self.v.t[1];
        let dz = self.u.t[2] - self.v.t[2];

        if dx < dy && dx < dz {
            0
        } else if dy < dz {
            1
        } else {
            2
        }
    }

    fn get_widest_dimension(&self) -> usize {
        let dx = self.u.t[0] - self.v.t[0];
        let dy = self.u.t[1] - self.v.t[1];
        let dz = self.u.t[2] - self.v.t[2];

        if dx > dy && dx > dz {
            0
        } else if dy > dz {
            1
        } else {
            2
        }
    }
}

impl BoundingBoxTrait for BoundingBox3d {
    fn partition(&self) -> Option<(Self, Self)> {
        if (self.u - self.v).length() < 1e-8 {
            return None;
        }

        let dim = self.get_widest_dimension();
        let midpoint = 0.5 * (self.v.t[dim] - self.u.t[dim]);

        let mut u_new = self.u;
        let mut v_new = self.v;
        u_new.t[dim] = midpoint;
        v_new.t[dim] = midpoint;

        Some((
            BoundingBox3d {
                u: self.u,
                v: v_new,
            },
            BoundingBox3d {
                u: u_new,
                v: self.v,
            },
        ))
    }

    fn extend(&self) -> Option<(Self, Self)> {
        if (self.u - self.v).length() > 1e8 {
            return None;
        }

        let dim = self.get_most_narrow_dimension();

        let right = self.v.t[dim];
        let left = self.u.t[dim];
        let dif = right - left;
        let left_clamped = if left < 0.0 { left } else { 0.0 };
        let right_clamped = if right > 0.0 { right } else { 0.0 };
        let mut template = self.clone();

        if -left_clamped > right_clamped {
            template.u.t[dim] -= dif;
            template.v.t[dim] -= dif;
            let parent = Self{u : template.u.clone(), v: self.v.clone()};
            Some((template, parent))
        } else {
            template.u.t[dim] += dif;
            template.v.t[dim] += dif;
            let parent = Self{u : self.u.clone(), v: template.v.clone()};
            Some((template, parent))
        }
    }

    fn is_contained(&self, other: &Self) -> bool {
        other.u.t[0] <= self.u.t[0]
            && self.v.t[0] <= other.v.t[0]
            && other.u.t[1] <= self.u.t[1]
            && self.v.t[1] <= other.v.t[1]
            && other.u.t[2] <= self.u.t[2]
            && self.v.t[2] <= other.v.t[2]
    }

    fn intersects(&self, other: &Self) -> bool {
        self.intersects_in_dimentions(other, 0)
            && self.intersects_in_dimentions(other, 1)
            && self.intersects_in_dimentions(other, 2)
    }

    fn is_sub_scale(&self, other: &Self) -> bool {
        (self.u - self.v).length() <= (other.u - other.v).length()
    }
}

impl HittableBoundingBoxTrait for BoundingBox3d {
    fn hit(&self, ray: &Ray, r: (f64, f64)) -> bool {
        let ray = ConstrainedRay{ray: ray.clone(), range : r};
        self.intersects(&ray).is_some()
    }
}

pub trait KdTreeContent<BoundingBox: BoundingBoxTrait>: Sized {
    fn get_bounding_box(&self) -> BoundingBox;
}

struct KdNode<BoundingBox: HittableBoundingBoxTrait, Content: KdTreeContent<BoundingBox>> {
    enclosure: BoundingBox,
    content: Vec<Arc<Content>>,
    children: Vec<Self>,
}

pub struct KdTree<BoundingBox: HittableBoundingBoxTrait, Content: KdTreeContent<BoundingBox>> {
    root: Option<KdNode<BoundingBox, Content>>,
}

struct ArcWrapper<T: ?Sized>(Arc<T>);

impl<T: ?Sized> std::cmp::PartialEq<ArcWrapper<T>> for ArcWrapper<T> {
    fn eq(&self, other: &Self) -> bool {
        Arc::as_ptr(&self.0) == Arc::as_ptr(&other.0)
    }
}

impl<T: ?Sized> std::cmp::Eq for ArcWrapper<T> {}

impl<T: ?Sized> Hash for ArcWrapper<T> {
    fn hash<H: Hasher>(&self, state: &mut H) {
        Arc::as_ptr(&self.0).hash(state);
    }
}

impl<BoundingBox: HittableBoundingBoxTrait, Content: KdTreeContent<BoundingBox>>
    KdTree<BoundingBox, Content>
{
    fn new_node(enclosure: BoundingBox) -> KdNode<BoundingBox, Content> {
        KdNode::<BoundingBox, Content> {
            enclosure,
            content: vec![],
            children: vec![],
        }
    }

    fn add_to_node(
        node: &mut KdNode<BoundingBox, Content>,
        content: &Arc<Content>,
        content_enclosure: &BoundingBox,
    ) {
        if content_enclosure.is_sub_scale(&node.enclosure) && node.children.is_empty() {
            if let Some((left, right)) = node.enclosure.partition() {
                let mut left = Self::new_node(left);
                let mut right = Self::new_node(right);

                for content in &node.content {
                    let enclosure = content.get_bounding_box();
                    if enclosure.intersects(&left.enclosure) {
                        left.content.push(content.clone());
                    }
                    if enclosure.intersects(&right.enclosure) {
                        right.content.push(content.clone());
                    }
                }
                node.content.clear();
                node.children.push(left);
                node.children.push(right);
            }
        }

        if node.children.is_empty() {
            node.content.push(content.clone());
        } else {
            for child in &mut node.children {
                Self::add_to_node(child, content, content_enclosure);
            }
        }
    }

    pub fn new() -> Self {
        Self {
            root: Some(Self::new_node(BoundingBox::default())),
        }
    }

    pub fn add(&mut self, content: &Arc<Content>) {
        let content_enclosure = content.get_bounding_box();
        let content = content.clone();

        while !content_enclosure.is_contained(&self.root.as_ref().unwrap().enclosure) {
            let (other, parent) = self.root.as_ref().unwrap().enclosure.extend().unwrap();

            let other_node = Self::new_node(other);
            let current_root = self.root.take().unwrap();
            let mut new_root = Self::new_node(parent);
            new_root.children.push(other_node);
            new_root.children.push(current_root);

            self.root = Some(new_root);
        }

        Self::add_to_node(
            &mut self.root.as_mut().unwrap(),
            &content,
            &content_enclosure,
        );
    }

    fn collect_ray_hit_candidates(
        node: &KdNode<BoundingBox, Content>,
        ray: &Ray,
        limits: (f64, f64),
        result: &mut HashSet<ArcWrapper<Content>>,
    ) {
        if !node.enclosure.hit(ray, limits) {
            return;
        }

        for content in &node.content {
            result.insert(ArcWrapper(content.clone()));
        }

        for children in &node.children {
            Self::collect_ray_hit_candidates(children, ray, limits, result);
        }
    }

    pub fn get_ray_hit_candidates(&self, ray: &Ray, limits: (f64, f64)) -> Vec<Arc<Content>> {
        let mut result = HashSet::<ArcWrapper<Content>>::new();
        Self::collect_ray_hit_candidates(self.root.as_ref().unwrap(), ray, limits, &mut result);
        result.into_iter().map(|x| x.0).collect::<Vec<_>>()
    }
}
