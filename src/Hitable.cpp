//
// Created by LamWS on 2021/5/14.
//

#include "Hitable.h"


using namespace Eigen;

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
