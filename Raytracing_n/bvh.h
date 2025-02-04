#pragma once
#ifndef BVHH
#define BVHH
#include "hitable.h"
#include "mathf.h"
#include "aabb.h"
#include "hitable_list.h"

class bvh_node :public hitable {
public:
	bvh_node() {}
	bvh_node(hitable **l, int n, float time0, float time1);
	virtual bool hit(const ray& r, float tmin, float tmax, hit_record& rec, bool is_medium = false) const;
	virtual bool bounding_box(float t0, float t1, aabb& box) const;
	hitable *left;
	hitable *right;
	hitable_list *list;
	aabb box;
};

int box_x_compare(const void * a, const void * b) {
	aabb box_left, box_right;
	hitable *ah = *(hitable**)a;
	hitable *bh = *(hitable**)b;
	if (!ah->bounding_box(0, 0, box_left) || !bh->bounding_box(0, 0, box_right))
		std::cerr << "no bounding box in bvh_node constructor\n";
	if (box_left.Min().x() - box_right.Min().x() < 0.0)
		return -1;
	else
		return 1;
}

int box_y_compare(const void * a, const void * b) {
	aabb box_left, box_right;
	hitable *ah = *(hitable**)a;
	hitable *bh = *(hitable**)b;
	if (!ah->bounding_box(0, 0, box_left) || !bh->bounding_box(0, 0, box_right))
		std::cerr << "no bounding box in bvh_node constructor\n";
	if (box_left.Min().y() - box_right.Min().y() < 0.0)
		return -1;
	else
		return 1;
}

int box_z_compare(const void * a, const void * b) {
	aabb box_left, box_right;
	hitable *ah = *(hitable**)a;
	hitable *bh = *(hitable**)b;
	if (!ah->bounding_box(0, 0, box_left) || !bh->bounding_box(0, 0, box_right))
		std::cerr << "no bounding box in bvh_node constructor\n";
	if (box_left.Min().z() - box_right.Min().z() < 0.0)
		return -1;
	else
		return 1;
}


bool bvh_node::bounding_box(float t0, float t1, aabb& b)const
{
	b = box;
	return true;
}

bool bvh_node::hit(const ray& r, float t_min, float t_max, hit_record& rec, bool is_medium) const
{
	if (box.hit(r, t_min, t_max))
	{
		hit_record left_rec, right_rec;
		bool hit_left = left->hit(r, t_min, t_max, left_rec, is_medium);
		bool hit_right = right->hit(r, t_min, t_max, right_rec, is_medium);
		if (hit_left && hit_right)
		{
			if (left_rec.t < right_rec.t)
				rec = left_rec;
			else
				rec = right_rec;
			return true;
		}
		else if (hit_left)
		{
			rec = left_rec;
			return true;
		}
		else if (hit_right)
		{
			rec = right_rec;
			return true;
		}
		else
			return false;
	}
	else return false;
}


bvh_node::bvh_node(hitable **l, int n, float time0, float time1) {
	int axis = int(3 * drand48());
	if (axis == 0)
		qsort(l, n, sizeof(hitable *), box_x_compare);
	else if (axis == 1)
		qsort(l, n, sizeof(hitable *), box_y_compare);
	else
		qsort(l, n, sizeof(hitable *), box_z_compare);
	if (n == 1) {
		left = right = l[0];
	}
	else if (n == 2) {
		left = l[0];
		right = l[1];
	}
	else {
		left = new bvh_node(l, n / 2, time0, time1);
		right = new bvh_node(l + n / 2, n - n / 2, time0, time1);
	}
	aabb box_left, box_right;
	if (!left->bounding_box(time0, time1, box_left) || !right->bounding_box(time0, time1, box_right))
		std::cerr << "no bounding box in bvh_node constructor\n";
	box = surrounding_box(box_left, box_right);
}

#endif // !BVHH

