//
// Created by LamWS on 2021/5/14.
//

#include "Hitable.h"


using namespace Eigen;

Aabb surrounding_box(const Aabb &box0, const Aabb &box1) {
    Vector3d small(fmin(box0.get_min().x(), box1.get_min().x()),
                   fmin(box0.get_min().y(), box1.get_min().y()),
                   fmin(box0.get_min().z(), box1.get_min().z()));
    Vector3d large(fmax(box0.get_max().x(), box1.get_max().x()),
                   fmax(box0.get_max().y(), box1.get_max().y()),
                   fmax(box0.get_max().z(), box1.get_max().z()));
    return Aabb(small, large);
}

bool Sphere::hit(const Ray &r, double t_min, double t_max, hit_record &rec) const {
    rec.material = material;
    Vector3d oc = r.origin() - center;
    double a = r.direction().dot(r.direction());
    double b = 2 * r.direction().dot(oc);
    double c = oc.dot(oc) - radius * radius;
    auto d = b * b - 4 * a * c;
    if (d > 0) {
        double tmp = (-b - sqrt(d)) / a / 2.0;
        if (tmp < t_max && tmp > t_min) {
            rec.t = tmp;
            rec.p = r.point_at_parameter(tmp);
            rec.normal = (rec.p - center) / radius;
            return true;
        }
        tmp = (-b + sqrt(d)) / a / 2.0;
        if (tmp < t_max && tmp > t_min) {
            rec.t = tmp;
            rec.p = r.point_at_parameter(tmp);
            rec.normal = (rec.p - center) / radius;
            return true;
        }
    }
    return false;
}

bool Sphere::bounding_box(double t0, double t1, Aabb &box) const {
    box = Aabb(center - Vector3d(radius, radius, radius), center + Vector3d(radius, radius, radius));
    return true;
}

Hitable_list::Hitable_list(std::vector<Hitable *> l, int n) : list(std::move(l)), list_size(n) {

}

bool Hitable_list::hit(const Ray &r, double t_min, double t_max, hit_record &rec) const {
    hit_record tmp;
    bool hit_anything = false;
    double closest = t_max;
    for (auto hitable:list) {
        if (hitable->hit(r, t_min, closest, tmp)) {
            hit_anything = true;
            closest = tmp.t;
            rec = tmp;
        }
    }
    return hit_anything;
}

bool Hitable_list::bounding_box(double t0, double t1, Aabb &box) const {
    if (list.empty())return false;
    Aabb temp_box;
    bool first_hit = list[0]->bounding_box(t0, t1, temp_box);
    if (!first_hit)return false;
    else box = temp_box;
    for (int i = 1; i < list.size(); i++) {
        if (list[i]->bounding_box(t0, t1, temp_box))
            box = surrounding_box(box, temp_box);
        else
            return false;
    }
    return true;
}

bool MovingSphere::hit(const Ray &r, double t_min, double t_max, hit_record &rec) const {
    rec.material = material;
    Vector3d oc = r.origin() - center(r.time());
    double a = r.direction().dot(r.direction());
    double b = 2 * r.direction().dot(oc);
    double c = oc.dot(oc) - radius * radius;
    auto d = b * b - 4 * a * c;
    if (d > 0) {
        double tmp = (-b - sqrt(d)) / a / 2.0;
        if (tmp < t_max && tmp > t_min) {
            rec.t = tmp;
            rec.p = r.point_at_parameter(tmp);
            rec.normal = (rec.p - center(r.time())) / radius;
            return true;
        }
        tmp = (-b + sqrt(d)) / a / 2.0;
        if (tmp < t_max && tmp > t_min) {
            rec.t = tmp;
            rec.p = r.point_at_parameter(tmp);
            rec.normal = (rec.p - center(r.time())) / radius;
            return true;
        }
    }
    return false;
}

Eigen::Vector3d MovingSphere::center(double t) const {
    return center0 + ((t - time0) / (time1 - time0)) * (center1 - center0);
}


bool MovingSphere::bounding_box(double t0, double t1, Aabb &box) const {
    Aabb box0(center(t0) - Vector3d(radius, radius, radius), center(t0) + Vector3d(radius, radius, radius));
    Aabb box1(center(t1) - Vector3d(radius, radius, radius), center(t1) + Vector3d(radius, radius, radius));
    box = surrounding_box(box0, box1);
    return true;
}


BVHNode::BVHNode(std::vector<Hitable *> list, double t0, double t1) {
    int axis = rand() % 3;

    auto x_cmp = [](Hitable *a, Hitable *b) {
        Aabb box_left, box_right;
        if (!a->bounding_box(0, 0, box_left) || !b->bounding_box(0, 0, box_right)) {

        }
        if (box_left.get_min().x() - box_right.get_min().x() < 0.0)return false;
        return true;
    };
    auto y_cmp = [](Hitable *a, Hitable *b) {
        Aabb box_left, box_right;
        if (!a->bounding_box(0, 0, box_left) || !b->bounding_box(0, 0, box_right)) {

        }
        if (box_left.get_min().y() - box_right.get_min().y() < 0.0)return false;
        return true;
    };
    auto z_cmp = [](Hitable *a, Hitable *b) {
        Aabb box_left, box_right;
        if (!a->bounding_box(0, 0, box_left) || !b->bounding_box(0, 0, box_right)) {

        }
        if (box_left.get_min().z() - box_right.get_min().z() < 0.0)return false;
        return true;
    };
    if (axis == 0) {
        std::sort(list.begin(), list.end(), x_cmp);
    } else if (axis == 1) {
        std::sort(list.begin(), list.end(), y_cmp);
    } else {
        std::sort(list.begin(), list.end(), z_cmp);
    }
    if (list.size() == 1) {
        left = right = list[0];
    } else if (list.size() == 2) {
        left = list[0];
        right = list[1];
    }else{
        left = new BVH
    }
    Aabb box_left, box_right;
    if (!left->bounding_box(t0, t1, box_left) || !right->bounding_box(t0, t1, box_right)) {

    }
    box = surrounding_box(box_right, box_right);
}

bool BVHNode::hit(const Ray &r, double t_min, double t_max, hit_record &rec) const {
    if (box.hit(r, t_min, t_max)) {
        hit_record left_rec, right_rec;
        bool hit_left = left->hit(r, t_min, t_max, left_rec);
        bool hit_right = right->hit(r, t_min, t_max, right_rec);
        if (hit_left && hit_right) {
            if (left_rec.t < right_rec.t)
                rec = left_rec;
            else
                rec = right_rec;
            return true;
        } else if (hit_left) {
            rec = left_rec;
            return true;
        } else if (hit_right) {
            rec = right_rec;
            return true;
        }
        return false;
    }
    return false;
}

bool BVHNode::bounding_box(double t0, double t1, Aabb &b) const {
    b = box;
    return true;
}
