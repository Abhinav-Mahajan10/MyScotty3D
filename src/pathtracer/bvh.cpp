
#include "bvh.h"
#include "aggregate.h"
#include "instance.h"
#include "tri_mesh.h"

#include <algorithm>
#include <array>
#include <limits>
#include <stack>

namespace PT {

struct BVHBuildData {
	BVHBuildData(size_t start, size_t range, size_t dst) : start(start), range(range), node(dst) {
	}
	size_t start; ///< start index into the primitive array
	size_t range; ///< range of index into the primitive array
	size_t node;  ///< address to update
};

struct SAHBucketData {
	BBox bb;          ///< bbox of all primitives
	size_t num_prims; ///< number of primitives in the bucket
};

template<typename Primitive>
void BVH<Primitive>::build(std::vector<Primitive>&& prims, size_t max_leaf_size)
{
	//A3T3 - build a bvh
	nodes.clear();
	primitives = std::move(prims);

	root_idx = 0;
	if (primitives.empty())
		return;

	root_idx = new_node();
	std::vector<BVHBuildData> stack;
	stack.emplace_back(0, primitives.size(), root_idx);

	while (!stack.empty())
	{
		BVHBuildData cur = stack.back();
		stack.pop_back();

		size_t start = cur.start;
		size_t range = cur.range;
		size_t end = start + range;

		BBox prim_bbox;
		BBox centroid_bbox;
		for (size_t i = start; i < end; i++)
		{
			BBox b = primitives[i].bbox();
			prim_bbox.enclose(b);
			centroid_bbox.enclose(b.center());
		}

		nodes[cur.node].bbox = prim_bbox;
		nodes[cur.node].start = start;
		nodes[cur.node].size = range;

		if (range <= max_leaf_size)
		{
			nodes[cur.node].l = nodes[cur.node].r = cur.node;
			continue;
		}

		constexpr int BUCKETS = 12;
		float best_cost = std::numeric_limits<float>::infinity();
		int best_axis = -1;
		int best_split_bucket = -1;

		for (int axis = 0; axis < 3; axis++)
		{
			float cmin = centroid_bbox.min[axis];
			float cmax = centroid_bbox.max[axis];
			if (cmax - cmin <= EPS_F)
				continue;

			std::array<SAHBucketData, BUCKETS> buckets;
			for (auto& b : buckets)
			{
				b.bb.reset();
				b.num_prims = 0;
			}

			for (size_t i = start; i < end; i++)
			{
				float c = primitives[i].bbox().center()[axis];
				float rel = (c - cmin) / (cmax - cmin);
				int b = std::min(BUCKETS - 1, std::max(0, static_cast<int>(rel * BUCKETS)));
				buckets[b].num_prims++;
				buckets[b].bb.enclose(primitives[i].bbox());
			}

			std::array<BBox, BUCKETS> left_bb, right_bb;
			std::array<size_t, BUCKETS> left_count{}, right_count{};
			BBox run_left, run_right;
			size_t cnt_left = 0, cnt_right = 0;
			for (int i = 0; i < BUCKETS; i++)
			{
				if (buckets[i].num_prims)
					run_left.enclose(buckets[i].bb);
				cnt_left += buckets[i].num_prims;
				left_bb[i] = run_left;
				left_count[i] = cnt_left;
			}
			for (int i = BUCKETS - 1; i >= 0; i--)
			{
				if (buckets[i].num_prims)
					run_right.enclose(buckets[i].bb);
				cnt_right += buckets[i].num_prims;
				right_bb[i] = run_right;
				right_count[i] = cnt_right;
			}

			float parent_area = prim_bbox.surface_area();
			if (parent_area <= EPS_F)
				continue;

			for (int i = 0; i < BUCKETS - 1; i++)
			{
				size_t lcount = left_count[i];
				size_t rcount = right_count[i + 1];
				if (lcount == 0 || rcount == 0)
					continue;
				float cost = (left_bb[i].surface_area() * float(lcount) +
				              right_bb[i + 1].surface_area() * float(rcount)) / parent_area;
				if (cost < best_cost)
				{
					best_cost = cost;
					best_axis = axis;
					best_split_bucket = i;
				}
			}
		}

		if (best_axis == -1)
		{
			nodes[cur.node].l = nodes[cur.node].r = cur.node;
			continue;
		}

		float cmin = centroid_bbox.min[best_axis];
		float cmax = centroid_bbox.max[best_axis];
		auto split_it = std::partition(primitives.begin() + static_cast<std::ptrdiff_t>(start),
		                               primitives.begin() + static_cast<std::ptrdiff_t>(end),
		                               [&](const Primitive& p) {
				float c = p.bbox().center()[best_axis];
				float rel = (c - cmin) / (cmax - cmin);
				int b = std::min(BUCKETS - 1, std::max(0, static_cast<int>(rel * BUCKETS)));
				return b <= best_split_bucket;
			});

		size_t mid = static_cast<size_t>(split_it - primitives.begin());
		if (mid == start || mid == end)
		{
			mid = start + range / 2;
			std::nth_element(primitives.begin() + static_cast<std::ptrdiff_t>(start),
			                 primitives.begin() + static_cast<std::ptrdiff_t>(mid),
			                 primitives.begin() + static_cast<std::ptrdiff_t>(end),
			                 [&](const Primitive& a, const Primitive& b) {
				return a.bbox().center()[best_axis] < b.bbox().center()[best_axis];
			});
		}

		size_t lidx = new_node();
		size_t ridx = new_node();
		nodes[cur.node].l = lidx;
		nodes[cur.node].r = ridx;

		stack.emplace_back(mid, end - mid, ridx);
		stack.emplace_back(start, mid - start, lidx);
	}
}

template<typename Primitive> Trace BVH<Primitive>::hit(const Ray& ray) const
{
	//A3T3 - traverse your BVH
	Trace ret;
	if (nodes.empty())
		return ret;

	float closest = ray.dist_bounds.y;
	std::vector<size_t> stack;
	stack.push_back(root_idx);

	while (!stack.empty())
	{
		size_t idx = stack.back();
		stack.pop_back();

		const Node& node = nodes[idx];
		Vec2 node_times(ray.dist_bounds.x, closest);
		if (!node.bbox.hit(ray, node_times))
			continue;

		if (node.is_leaf())
		{
			for (size_t i = node.start; i < node.start + node.size; i++)
			{
				Ray local(ray.point, ray.dir, Vec2(ray.dist_bounds.x, closest), ray.depth);
				Trace h = primitives[i].hit(local);
				if (h.hit && h.distance < closest)
				{
					closest = h.distance;
					ret = h;
				}
			}
			continue;
		}

		Vec2 ltimes(ray.dist_bounds.x, closest), rtimes(ray.dist_bounds.x, closest);
		bool lhit = nodes[node.l].bbox.hit(ray, ltimes);
		bool rhit = nodes[node.r].bbox.hit(ray, rtimes);

		if (lhit && rhit)
		{
			if (ltimes.x < rtimes.x)
			{
				stack.push_back(node.r);
				stack.push_back(node.l);
			}
			else
			{
				stack.push_back(node.l);
				stack.push_back(node.r);
			}
		}
		else if (lhit)
			stack.push_back(node.l);
		else if (rhit)
			stack.push_back(node.r);
	}

	return ret;
}

template<typename Primitive>
BVH<Primitive>::BVH(std::vector<Primitive>&& prims, size_t max_leaf_size) {
	build(std::move(prims), max_leaf_size);
}

template<typename Primitive> std::vector<Primitive> BVH<Primitive>::destructure() {
	nodes.clear();
	return std::move(primitives);
}

template<typename Primitive>
template<typename P>
typename std::enable_if<std::is_copy_assignable_v<P>, BVH<P>>::type BVH<Primitive>::copy() const {
	BVH<Primitive> ret;
	ret.nodes = nodes;
	ret.primitives = primitives;
	ret.root_idx = root_idx;
	return ret;
}

template<typename Primitive> Vec3 BVH<Primitive>::sample(RNG &rng, Vec3 from) const {
	if (primitives.empty()) return {};
	int32_t n = rng.integer(0, static_cast<int32_t>(primitives.size()));
	return primitives[n].sample(rng, from);
}

template<typename Primitive>
float BVH<Primitive>::pdf(Ray ray, const Mat4& T, const Mat4& iT) const {
	if (primitives.empty()) return 0.0f;
	float ret = 0.0f;
	for (auto& prim : primitives) ret += prim.pdf(ray, T, iT);
	return ret / primitives.size();
}

template<typename Primitive> void BVH<Primitive>::clear() {
	nodes.clear();
	primitives.clear();
}

template<typename Primitive> bool BVH<Primitive>::Node::is_leaf() const {
	// A node is a leaf if l == r, since all interior nodes must have distinct children
	return l == r;
}

template<typename Primitive>
size_t BVH<Primitive>::new_node(BBox box, size_t start, size_t size, size_t l, size_t r) {
	Node n;
	n.bbox = box;
	n.start = start;
	n.size = size;
	n.l = l;
	n.r = r;
	nodes.push_back(n);
	return nodes.size() - 1;
}
 
template<typename Primitive> BBox BVH<Primitive>::bbox() const {
	if(nodes.empty()) return BBox{Vec3{0.0f}, Vec3{0.0f}};
	return nodes[root_idx].bbox;
}

template<typename Primitive> size_t BVH<Primitive>::n_primitives() const {
	return primitives.size();
}

template<typename Primitive>
uint32_t BVH<Primitive>::visualize(GL::Lines& lines, GL::Lines& active, uint32_t level,
                                   const Mat4& trans) const {

	std::stack<std::pair<size_t, uint32_t>> tstack;
	tstack.push({root_idx, 0u});
	uint32_t max_level = 0u;

	if (nodes.empty()) return max_level;

	while (!tstack.empty()) {

		auto [idx, lvl] = tstack.top();
		max_level = std::max(max_level, lvl);
		const Node& node = nodes[idx];
		tstack.pop();

		Spectrum color = lvl == level ? Spectrum(1.0f, 0.0f, 0.0f) : Spectrum(1.0f);
		GL::Lines& add = lvl == level ? active : lines;

		BBox box = node.bbox;
		box.transform(trans);
		Vec3 min = box.min, max = box.max;

		auto edge = [&](Vec3 a, Vec3 b) { add.add(a, b, color); };

		edge(min, Vec3{max.x, min.y, min.z});
		edge(min, Vec3{min.x, max.y, min.z});
		edge(min, Vec3{min.x, min.y, max.z});
		edge(max, Vec3{min.x, max.y, max.z});
		edge(max, Vec3{max.x, min.y, max.z});
		edge(max, Vec3{max.x, max.y, min.z});
		edge(Vec3{min.x, max.y, min.z}, Vec3{max.x, max.y, min.z});
		edge(Vec3{min.x, max.y, min.z}, Vec3{min.x, max.y, max.z});
		edge(Vec3{min.x, min.y, max.z}, Vec3{max.x, min.y, max.z});
		edge(Vec3{min.x, min.y, max.z}, Vec3{min.x, max.y, max.z});
		edge(Vec3{max.x, min.y, min.z}, Vec3{max.x, max.y, min.z});
		edge(Vec3{max.x, min.y, min.z}, Vec3{max.x, min.y, max.z});

		if (!node.is_leaf()) {
			tstack.push({node.l, lvl + 1});
			tstack.push({node.r, lvl + 1});
		} else {
			for (size_t i = node.start; i < node.start + node.size; i++) {
				uint32_t c = primitives[i].visualize(lines, active, level - lvl, trans);
				max_level = std::max(c + lvl, max_level);
			}
		}
	}
	return max_level;
}

template class BVH<Triangle>;
template class BVH<Instance>;
template class BVH<Aggregate>;
template BVH<Triangle> BVH<Triangle>::copy<Triangle>() const;

} // namespace PT
