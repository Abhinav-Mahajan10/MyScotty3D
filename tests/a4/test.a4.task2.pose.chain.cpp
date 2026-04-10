#include "test.h"
#include "scene/skeleton.h"

Test test_a4_task2_pose_bind_chain("a4.task2.pose.bind.chain", []()
{
	Skeleton sk;
	sk.base = Vec3(0.0f, 1.0f, 0.0f);
	auto root = sk.add_bone(-1U, Vec3(2.0f, 0.0f, 0.0f));
	auto mid = sk.add_bone(root, Vec3(0.0f, 3.0f, 0.0f));
	auto tip_bone = sk.add_bone(mid, Vec3(0.0f, 0.0f, 1.0f));

	(void)tip_bone;

	std::vector<Mat4> B = sk.bind_pose();

	Vec3 root_o = B[root] * Vec3(0.0f, 0.0f, 0.0f);
	Vec3 mid_o = B[mid] * Vec3(0.0f, 0.0f, 0.0f);
	Vec3 tip = B[tip_bone] * Vec3(0.0f, 0.0f, 0.0f);

	if (Test::differs(root_o, Vec3(0.0f, 1.0f, 0.0f)))
	{
		throw Test::error("Root bone origin mismatch");
	}
	if (Test::differs(mid_o, Vec3(2.0f, 1.0f, 0.0f)))
	{
		throw Test::error("Mid bone origin mismatch");
	}
	if (Test::differs(tip, Vec3(2.0f, 4.0f, 0.0f)))
	{
		throw Test::error("Tip bone origin mismatch");
	}
});

Test test_a4_task2_base_offset_current("a4.task2.pose.current.base_offset", []()
{
	Skeleton sk;
	sk.base = Vec3(0.0f, 0.0f, 0.0f);
	sk.base_offset = Vec3(0.0f, 0.0f, 5.0f);
	auto b = sk.add_bone(-1U, Vec3(1.0f, 0.0f, 0.0f));
	std::vector<Mat4> P = sk.current_pose();
	Vec3 o = P[b] * Vec3(0.0f, 0.0f, 0.0f);
	if (Test::differs(o, Vec3(0.0f, 0.0f, 5.0f)))
	{
		throw Test::error("base_offset should shift root bone origin in current_pose");
	}
});

Test test_a4_task2_bind_multi_child("a4.task2.pose.bind.multi_child", []()
{
	Skeleton sk;
	sk.base = Vec3(0.0f, 0.0f, 0.0f);
	auto root  = sk.add_bone(-1U,  Vec3(1.0f, 0.0f, 0.0f));
	auto childA = sk.add_bone(root, Vec3(0.0f, 1.0f, 0.0f));
	auto childB = sk.add_bone(root, Vec3(0.0f, 0.0f, 1.0f));
	(void)childA; (void)childB;

	std::vector<Mat4> B = sk.bind_pose();

	Vec3 oA = B[childA] * Vec3(0.0f, 0.0f, 0.0f);
	Vec3 oB = B[childB] * Vec3(0.0f, 0.0f, 0.0f);
	if (Test::differs(oA, Vec3(1.0f, 0.0f, 0.0f)))
	{
		throw Test::error("childA origin should be at parent tip (1,0,0)");
	}
	if (Test::differs(oB, Vec3(1.0f, 0.0f, 0.0f)))
	{
		throw Test::error("childB origin should be at parent tip (1,0,0)");
	}
});

Test test_a4_task2_ik_base_offset_unchanged("a4.task2.step_ik.base_offset_unchanged", []()
{
	Skeleton sk;
	sk.base        = Vec3(0.0f, 0.0f, 0.0f);
	sk.base_offset = Vec3(1.0f, 2.0f, 3.0f);
	auto bone   = sk.add_bone(-1U, Vec3(0.0f, 1.0f, 0.0f));
	auto handle = sk.add_handle(bone, Vec3(0.0f, 0.5f, 0.5f));
	sk.handles[handle].enabled = true;

	sk.solve_ik(200);

	if (Test::differs(sk.base_offset, Vec3(1.0f, 2.0f, 3.0f)))
	{
		throw Test::error("solve_ik must not modify base_offset");
	}
});
