#include <algorithm>
#include <unordered_set>
#include "skeleton.h"
#include "test.h"
#include <iostream>

void Skeleton::Bone::compute_rotation_axes(Vec3 *x_, Vec3 *y_, Vec3 *z_) const {
	assert(x_ && y_ && z_);
	auto &x = *x_;
	auto &y = *y_;
	auto &z = *z_;

	//y axis points in the direction of extent:
	y = extent.unit();
	//if extent is too short to normalize nicely, point along the skeleton's 'y' axis:
	if (!y.valid()) {
		y = Vec3{0.0f, 1.0f, 0.0f};
	}

	//x gets skeleton's 'x' axis projected to be orthogonal to 'y':
	x = Vec3{1.0f, 0.0f, 0.0f};
	x = (x - dot(x,y) * y).unit();
	if (!x.valid()) {
		//if y perfectly aligns with skeleton's 'x' axis, x, gets skeleton's z axis:
		x = Vec3{0.0f, 0.0f, 1.0f};
		x = (x - dot(x,y) * y).unit(); //(this should do nothing)
	}

	//z computed from x,y:
	z = cross(x,y);

	//x,z rotated by roll:
	float cr = std::cos(roll / 180.0f * PI_F);
	float sr = std::sin(roll / 180.0f * PI_F);
	// x = cr * x + sr * -z;
	// z = cross(x,y);
	std::tie(x, z) = std::make_pair(cr * x + sr * -z, cr * z + sr * x);
}

std::vector< Mat4 > Skeleton::bind_pose() const {
	//A4T2a: bone-to-skeleton transformations in the bind pose
	//(the bind pose does not rotate by Bone::pose)

	std::vector< Mat4 > bind;
	bind.reserve(bones.size());

	//NOTE: bones is guaranteed to be ordered such that parents appear before child bones.

	for (uint32_t bi = 0; bi < bones.size(); ++bi)
	{
		auto const &bone = bones[bi];
		(void)bone; //avoid complaints about unused bone
		//placeholder -- your code should actually compute the correct transform:
		if (bone.parent == -1U)
		{
			bind.push_back(Mat4::translate(base));
		}
		else
		{
			bind.push_back(bind[bone.parent] * Mat4::translate(bones[bone.parent].extent));
		}
	}

	assert(bind.size() == bones.size()); //should have a transform for every bone.
	return bind;
}

std::vector< Mat4 > Skeleton::current_pose() const {
    //A4T2a: bone-to-skeleton transformations in the current pose

	//Similar to bind_pose(), but takes rotation from Bone::pose into account.
	// (and translation from Skeleton::base_offset!)

	//You'll probably want to write a loop similar to bind_pose().

	//Useful functions:
	//Bone::compute_rotation_axes() will tell you what axes (in local bone space) Bone::pose should rotate around.
	//Mat4::angle_axis(angle, axis) will produce a matrix that rotates angle (in degrees) around a given axis.

	std::vector< Mat4 > pose;
	pose.reserve(bones.size());

	for (uint32_t bi = 0; bi < bones.size(); ++bi)
	{
		auto const &bone = bones[bi];
		Vec3 ax, ay, az;
		bone.compute_rotation_axes(&ax, &ay, &az);
		Mat4 R = Mat4::angle_axis(bone.pose.z, az) * Mat4::angle_axis(bone.pose.y, ay) *
		         Mat4::angle_axis(bone.pose.x, ax);
		if (bone.parent == -1U)
		{
			pose.push_back(Mat4::translate(base + base_offset) * R);
		}
		else
		{
			pose.push_back(pose[bone.parent] * Mat4::translate(bones[bone.parent].extent) * R);
		}
	}

	return pose;
}

std::vector< Vec3 > Skeleton::gradient_in_current_pose() const {
    //A4T2b: IK gradient

    // Computes the gradient (partial derivative) of IK energy relative to each bone's Bone::pose, in the current pose.

	//The IK energy is the sum over all *enabled* handles of the squared distance from the tip of Handle::bone to Handle::target
	std::vector< Vec3 > gradient(bones.size(), Vec3{0.0f, 0.0f, 0.0f});

	//TODO: loop over handles and over bones in the chain leading to the handle, accumulating gradient contributions.
	//remember bone.compute_rotation_axes() -- should be useful here, too!

	std::vector< Mat4 > pose = current_pose();

	for (Handle const &handle : handles)
	{
		if (!handle.enabled)
		{
			continue;
		}
		BoneIndex hi = handle.bone;
		Vec3 tip = pose[hi] * bones[hi].extent;
		Vec3 err = tip - handle.target;

		for (BoneIndex b = hi; b != -1U; b = bones[b].parent)
		{
			Bone const &bone = bones[b];
			Vec3 ax, ay, az;
			bone.compute_rotation_axes(&ax, &ay, &az);
			Mat4 R_x = Mat4::angle_axis(bone.pose.x, ax);
			Mat4 R_y = Mat4::angle_axis(bone.pose.y, ay);
			Mat4 R_z = Mat4::angle_axis(bone.pose.z, az);

			Mat4 M_joint;
			if (bone.parent == -1U)
			{
				M_joint = Mat4::translate(base + base_offset);
			}
			else
			{
				M_joint = pose[bone.parent] * Mat4::translate(bones[bone.parent].extent);
			}

			Vec3 r = M_joint * Vec3(0.0f, 0.0f, 0.0f);

			Mat4 M_rx = M_joint * R_z * R_y * R_x;
			Vec3 axis_rx = M_rx.rotate(Vec3(1.0f, 0.0f, 0.0f));
			Vec3 d_rx = cross(axis_rx, tip - r);

			Mat4 M_ry = M_joint * R_z * R_y;
			Vec3 axis_ry = M_ry.rotate(Vec3(0.0f, 1.0f, 0.0f));
			Vec3 d_ry = cross(axis_ry, tip - r);

			Mat4 M_rz = M_joint * R_z;
			Vec3 axis_rz = M_rz.rotate(Vec3(0.0f, 0.0f, 1.0f));
			Vec3 d_rz = cross(axis_rz, tip - r);

			gradient[b].x += dot(err, d_rx);
			gradient[b].y += dot(err, d_ry);
			gradient[b].z += dot(err, d_rz);
		}
	}

	assert(gradient.size() == bones.size());
	return gradient;
}

bool Skeleton::solve_ik(uint32_t steps) {
	//A4T2b - gradient descent
	//check which handles are enabled
	//run `steps` iterations

	//call gradient_in_current_pose() to compute d loss / d pose
	//add ...

	//if at a local minimum (e.g., gradient is near-zero), return 'true'.
	//if run through all steps, return `false`.

	const float grad_eps = 1.0e-4f;

	auto compute_cost = [&]() -> float
	{
		std::vector< Mat4 > pose = current_pose();
		float cost = 0.0f;
		for (Handle const &h : handles)
		{
			if (!h.enabled) continue;
			Vec3 tip = pose[h.bone] * bones[h.bone].extent;
			Vec3 err = tip - h.target;
			cost += 0.5f * dot(err, err);
		}
		return cost;
	};

	for (uint32_t s = 0; s < steps; ++s)
	{
		std::vector< Vec3 > grad = gradient_in_current_pose();
		float g2 = 0.0f;
		for (Vec3 const &g : grad)
		{
			g2 += dot(g, g);
		}
		if (g2 < grad_eps * grad_eps)
		{
			return true;
		}

		float f0 = compute_cost();

		std::vector< Vec3 > saved(bones.size());
		for (uint32_t bi = 0; bi < bones.size(); ++bi)
		{
			saved[bi] = bones[bi].pose;
		}

		float tau = 1.0f;
		bool accepted = false;
		for (int ls = 0; ls < 32; ++ls)
		{
			for (uint32_t bi = 0; bi < bones.size(); ++bi)
			{
				bones[bi].pose = saved[bi] - tau * grad[bi];
			}
			if (compute_cost() < f0)
			{
				accepted = true;
				break;
			}
			tau *= 0.5f;
		}

		if (!accepted)
		{
			for (uint32_t bi = 0; bi < bones.size(); ++bi)
			{
				bones[bi].pose = saved[bi];
			}
			return true;
		}
	}

	return false;
}

Vec3 Skeleton::closest_point_on_line_segment(Vec3 const &a, Vec3 const &b, Vec3 const &p)
{
	//A4T3: bone weight computation (closest point helper)
	//
	// Return the closest point to 'p' on the line segment from a to b
	//
	// Efficiency note: you can do this without any sqrt's! (no .unit() or .norm() is needed!)

	Vec3 ab = b - a;
	float denom = dot(ab, ab);
	if (denom <= EPS_F * EPS_F)
	{
		return a;
	}
	float t = dot(p - a, ab) / denom;
	t = std::clamp(t, 0.0f, 1.0f);
	return a + t * ab;
}

void Skeleton::assign_bone_weights(Halfedge_Mesh *mesh_) const
{
	assert(mesh_);
	auto &mesh = *mesh_;

	//A4T3: bone weight computation
	//
	// Visit every vertex and **set new values** in Vertex::bone_weights (don't append to old values)
	//
	// Be sure to use bone positions in the bind pose (not the current pose!)
	//
	// You should fill in the helper closest_point_on_line_segment() before working on this function

	std::vector< Mat4 > bind = bind_pose();

	for (auto vi = mesh.vertices.begin(); vi != mesh.vertices.end(); ++vi)
	{
		vi->bone_weights.clear();
		float sum = 0.0f;
		std::vector< Halfedge_Mesh::Vertex::Bone_Weight > raw;
		raw.reserve(bones.size());
		Vec3 const &pos = vi->position;

		for (uint32_t j = 0; j < bones.size(); ++j)
		{
			float r = bones[j].radius;
			if (r <= EPS_F)
			{
				continue;
			}
			Vec3 seg_a = bind[j] * Vec3(0.0f, 0.0f, 0.0f);
			Vec3 seg_b = bind[j] * bones[j].extent;
			Vec3 c = closest_point_on_line_segment(seg_a, seg_b, pos);
			float d = (pos - c).norm();
			float w_hat = (r - d) / r;
			if (w_hat > 0.0f)
			{
				raw.push_back(Halfedge_Mesh::Vertex::Bone_Weight{j, w_hat});
				sum += w_hat;
			}
		}

		if (sum <= EPS_F)
		{
			continue;
		}

		for (auto const &bw : raw)
		{
			vi->bone_weights.push_back(
				Halfedge_Mesh::Vertex::Bone_Weight{bw.bone, bw.weight / sum});
		}
	}
}

Indexed_Mesh Skeleton::skin(Halfedge_Mesh const &mesh, std::vector< Mat4 > const &bind,
                            std::vector< Mat4 > const &current)
{
	assert(bind.size() == current.size());

	//A4T3: linear blend skinning
	//
	// One approach you might take is to first compute the skinned positions (at every vertex) and normals (at every corner)
	// then generate faces in the style of Indexed_Mesh::from_halfedge_mesh
	//
	// ---- step 1: figure out skinned positions ---
	// (you will probably want to precompute some bind-to-current transformation matrices here)
	//
	// ---- step 2: transform into an indexed mesh ---
	// Hint: you should be able to use the code from Indexed_Mesh::from_halfedge_mesh (SplitEdges version) pretty much verbatim,
	// you'll just need to fill in the positions and normals.
	//
	// Skinned vertex: v' = sum_j w_ij P_j B_j^{-1} v

	uint32_t const n_bones = static_cast<uint32_t>(bind.size());
	std::vector< Mat4 > skin_M;
	skin_M.reserve(n_bones);
	for (uint32_t j = 0; j < n_bones; ++j)
	{
		skin_M.push_back(current[j] * bind[j].inverse());
	}

	std::unordered_map< Halfedge_Mesh::VertexCRef, Vec3 > skinned_positions;
	std::unordered_map< Halfedge_Mesh::HalfedgeCRef, Vec3 > skinned_normals;
	skinned_positions.reserve(mesh.vertices.size());
	skinned_normals.reserve(mesh.halfedges.size());

	for (auto vi = mesh.vertices.begin(); vi != mesh.vertices.end(); ++vi)
	{
		Vec3 pos;
		Mat4 T = Mat4::Zero;
		if (vi->bone_weights.empty())
		{
			pos = vi->position;
		}
		else
		{
			for (auto const &bw : vi->bone_weights)
			{
				T += skin_M[bw.bone] * bw.weight;
			}
			pos = T * vi->position;
		}
		skinned_positions.emplace(vi, pos);

		Mat4 Nxf = Mat4::I;
		if (!vi->bone_weights.empty())
		{
			Nxf = T.inverse().T();
		}

		auto h = vi->halfedge;
		do
		{
			if (!h->face->boundary)
			{
				Vec3 n = Nxf.rotate(h->corner_normal).unit();
				skinned_normals.emplace(h, n);
			}
			h = h->twin->next;
		} while (h != vi->halfedge);
	}

	std::vector< Indexed_Mesh::Vert > verts;
	std::vector< Indexed_Mesh::Index > idxs;

	for (Halfedge_Mesh::FaceCRef f = mesh.faces.begin(); f != mesh.faces.end(); f++)
	{
		if (f->boundary)
			continue;

		uint32_t corners_begin = static_cast<uint32_t>(verts.size());
		Halfedge_Mesh::HalfedgeCRef h = f->halfedge;
		do
		{
			Indexed_Mesh::Vert vert;
			auto pit = skinned_positions.find(h->vertex);
			auto nit = skinned_normals.find(h);
			assert(pit != skinned_positions.end());
			assert(nit != skinned_normals.end());
			vert.pos = pit->second;
			vert.norm = nit->second;
			vert.uv = h->corner_uv;
			vert.id = f->id;
			verts.emplace_back(vert);
			h = h->next;
		} while (h != f->halfedge);
		uint32_t const corners_end = static_cast<uint32_t>(verts.size());
		for (size_t i = corners_begin + 1; i + 1 < corners_end; i++)
		{
			idxs.emplace_back(corners_begin);
			idxs.emplace_back(static_cast<uint32_t>(i));
			idxs.emplace_back(static_cast<uint32_t>(i + 1));
		}
	}

	return Indexed_Mesh(std::move(verts), std::move(idxs));
}

void Skeleton::for_bones(const std::function<void(Bone&)>& f) {
	for (auto& bone : bones) {
		f(bone);
	}
}


void Skeleton::erase_bone(BoneIndex bone) {
	assert(bone < bones.size());
	//update indices in bones:
	for (uint32_t b = 0; b < bones.size(); ++b) {
		if (bones[b].parent == -1U) continue;
		if (bones[b].parent == bone) {
			assert(b > bone); //topological sort!
			//keep bone tips in the same place when deleting parent bone:
			bones[b].extent += bones[bone].extent;
			bones[b].parent = bones[bone].parent;
		} else if (bones[b].parent > bone) {
			assert(b > bones[b].parent); //topological sort!
			bones[b].parent -= 1;
		}
	}
	// erase the bone
	bones.erase(bones.begin() + bone);
	//update indices in handles (and erase any handles on this bone):
	for (uint32_t h = 0; h < handles.size(); /* later */) {
		if (handles[h].bone == bone) {
			erase_handle(h);
		} else if (handles[h].bone > bone) {
			handles[h].bone -= 1;
			++h;
		} else {
			++h;
		}
	}
}

void Skeleton::erase_handle(HandleIndex handle) {
	assert(handle < handles.size());

	//nothing internally refers to handles by index so can just delete:
	handles.erase(handles.begin() + handle);
}


Skeleton::BoneIndex Skeleton::add_bone(BoneIndex parent, Vec3 extent) {
	assert(parent == -1U || parent < bones.size());
	Bone bone;
	bone.extent = extent;
	bone.parent = parent;
	//all other parameters left as default.

	//slightly unfortunate hack:
	//(to ensure increasing IDs within an editing session, but reset on load)
	std::unordered_set< uint32_t > used;
	for (auto const &b : bones) {
		used.emplace(b.channel_id);
	}
	while (used.count(next_bone_channel_id)) ++next_bone_channel_id;
	bone.channel_id = next_bone_channel_id++;

	//all other parameters left as default.

	BoneIndex index = BoneIndex(bones.size());
	bones.emplace_back(bone);

	return index;
}

Skeleton::HandleIndex Skeleton::add_handle(BoneIndex bone, Vec3 target) {
	assert(bone < bones.size());
	Handle handle;
	handle.bone = bone;
	handle.target = target;
	//all other parameters left as default.

	//slightly unfortunate hack:
	//(to ensure increasing IDs within an editing session, but reset on load)
	std::unordered_set< uint32_t > used;
	for (auto const &h : handles) {
		used.emplace(h.channel_id);
	}
	while (used.count(next_handle_channel_id)) ++next_handle_channel_id;
	handle.channel_id = next_handle_channel_id++;

	HandleIndex index = HandleIndex(handles.size());
	handles.emplace_back(handle);

	return index;
}


Skeleton Skeleton::copy() {
	//turns out that there aren't any fancy pointer data structures to fix up here.
	return *this;
}

void Skeleton::make_valid() {
	for (uint32_t b = 0; b < bones.size(); ++b) {
		if (!(bones[b].parent == -1U || bones[b].parent < b)) {
			warn("bones[%u].parent is %u, which is not < %u; setting to -1.", b, bones[b].parent, b);
			bones[b].parent = -1U;
		}
	}
	if (bones.empty() && !handles.empty()) {
		warn("Have %u handles but no bones. Deleting handles.", uint32_t(handles.size()));
		handles.clear();
	}
	for (uint32_t h = 0; h < handles.size(); ++h) {
		if (handles[h].bone >= HandleIndex(bones.size())) {
			warn("handles[%u].bone is %u, which is not < bones.size(); setting to 0.", h, handles[h].bone);
			handles[h].bone = 0;
		}
	}
}

//-------------------------------------------------

Indexed_Mesh Skinned_Mesh::bind_mesh() const {
	return Indexed_Mesh::from_halfedge_mesh(mesh, Indexed_Mesh::SplitEdges);
}

Indexed_Mesh Skinned_Mesh::posed_mesh() const {
	return Skeleton::skin(mesh, skeleton.bind_pose(), skeleton.current_pose());
}

Skinned_Mesh Skinned_Mesh::copy() {
	return Skinned_Mesh{mesh.copy(), skeleton.copy()};
}
