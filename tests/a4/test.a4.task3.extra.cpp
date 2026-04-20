#include "test.h"
#include "geometry/util.h"
#include "scene/skeleton.h"

Test test_a4_task3_closest_point_off_axis(
	"a4.task3.closest_point_on_line_segment.off_axis", []()
	{
		Vec3 a = Vec3(0.0f, 0.0f, 0.0f);
		Vec3 b = Vec3(0.0f, 2.0f, 0.0f);
		Vec3 p = Vec3(3.0f, 1.0f, 0.0f);
		Vec3 expected = Vec3(0.0f, 1.0f, 0.0f);
		Vec3 actual = Skeleton::closest_point_on_line_segment(a, b, p);
		if (Test::differs(expected, actual))
		{
			throw Test::error("Off-axis point should project to segment midpoint.");
		}
	});

Test test_a4_task3_closest_point_degenerate(
	"a4.task3.closest_point_on_line_segment.degenerate", []()
	{
		Vec3 a = Vec3(5.0f, 3.0f, 1.0f);
		Vec3 b = Vec3(5.0f, 3.0f, 1.0f);
		Vec3 p = Vec3(0.0f, 0.0f, 0.0f);
		Vec3 actual = Skeleton::closest_point_on_line_segment(a, b, p);
		if (Test::differs(a, actual))
		{
			throw Test::error("Degenerate segment should return the single endpoint.");
		}
	});

Test test_a4_task3_assign_bone_weights_outside_radius(
	"a4.task3.assign_bone_weights.outside_radius", []()
	{
		Indexed_Mesh cube = Util::cube_mesh(0.1f);
		Halfedge_Mesh mesh = Halfedge_Mesh::from_indexed_mesh(cube);

		Skeleton sk;
		sk.base = Vec3(100.0f, 0.0f, 0.0f);
		auto b = sk.add_bone(-1U, Vec3(1.0f, 0.0f, 0.0f));
		sk.bones[b].radius = 0.5f;

		sk.assign_bone_weights(&mesh);

		for (auto vi = mesh.vertices.begin(); vi != mesh.vertices.end(); ++vi)
		{
			if (!vi->bone_weights.empty())
			{
				throw Test::error("Vertex outside all bone radii should have no weights.");
			}
		}
	});

Test test_a4_task3_assign_bone_weights_equidistant(
	"a4.task3.assign_bone_weights.equidistant_two_bones", []()
	{
		Skeleton sk;
		sk.base = Vec3(0.0f, 0.0f, 0.0f);
		auto b1 = sk.add_bone(-1U, Vec3(0.0f, 1.0f, 0.0f));
		sk.bones[b1].radius = 2.0f;

		sk.base = Vec3(2.0f, 0.0f, 0.0f);
		auto b2 = sk.add_bone(-1U, Vec3(0.0f, 1.0f, 0.0f));
		sk.bones[b2].radius = 2.0f;

		Halfedge_Mesh mesh;
		auto v = mesh.emplace_vertex();
		v->position = Vec3(1.0f, 0.5f, 0.0f);

		sk.assign_bone_weights(&mesh);

		if (v->bone_weights.size() != 2)
		{
			throw Test::error("Equidistant vertex should have weights for both bones.");
		}
		if (std::abs(v->bone_weights[0].weight - 0.5f) > 1e-5f ||
		    std::abs(v->bone_weights[1].weight - 0.5f) > 1e-5f)
		{
			throw Test::error("Equidistant vertex should get equal weight (0.5) from each bone.");
		}
	});

Test test_a4_task3_skin_identity_pose(
	"a4.task3.skin.identity_pose", []()
	{
		Indexed_Mesh cyl = Util::cyl_mesh(0.5f, 2);
		Halfedge_Mesh mesh = Halfedge_Mesh::from_indexed_mesh(cyl);

		Skeleton sk;
		sk.base = Vec3(0.0f, 0.0f, 0.0f);
		auto b = sk.add_bone(-1U, Vec3(0.0f, 2.0f, 0.0f));
		sk.bones[b].radius = 1.0f;

		auto bind = sk.bind_pose();
		auto cur = sk.current_pose();
		sk.assign_bone_weights(&mesh);
		auto skinned = sk.skin(mesh, bind, cur);

		std::vector<Vec3> orig;
		for (auto vi = mesh.vertices.begin(); vi != mesh.vertices.end(); ++vi)
		{
			orig.push_back(vi->position);
		}

		for (auto const &sv : skinned.vertices())
		{
			bool found = false;
			for (auto const &op : orig)
			{
				if (!Test::differs(sv.pos, op))
				{
					found = true;
					break;
				}
			}
			if (!found)
			{
				throw Test::error("skin() with identity pose should leave all vertex positions unchanged.");
			}
		}
	});
