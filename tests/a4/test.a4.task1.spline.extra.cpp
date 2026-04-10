#include "test.h"
#include "geometry/spline.h"

Test test_a4_task1_spline_empty("a4.task1.spline.extra.empty", []()
{
	Spline<Vec3> s;
	if (Test::differs(s.at(0.0f), Vec3(0.0f, 0.0f, 0.0f)))
	{
		throw Test::error("Empty spline should return default-constructed value");
	}
});

Test test_a4_task1_spline_single_knot("a4.task1.spline.extra.single_knot", []()
{
	Spline<Vec3> s;
	s.set(3.0f, Vec3(1.0f, 2.0f, 3.0f));
	Vec3 v = s.at(3.0f);
	if (Test::differs(v, Vec3(1.0f, 2.0f, 3.0f)))
	{
		throw Test::error("Spline should hit the single knot at its time");
	}
	v = s.at(-100.0f);
	if (Test::differs(v, Vec3(1.0f, 2.0f, 3.0f)))
	{
		throw Test::error("Single-knot spline should be constant in time");
	}
});

Test test_a4_task1_spline_clamp_ends("a4.task1.spline.extra.clamp_ends", []()
{
	Spline<float> s;
	s.set(1.0f, 10.0f);
	s.set(4.0f, 40.0f);
	s.set(9.0f, 90.0f);
	if (Test::differs(s.at(0.5f), 10.0f))
	{
		throw Test::error("Should clamp to first knot before its time");
	}
	if (Test::differs(s.at(100.0f), 90.0f))
	{
		throw Test::error("Should clamp to last knot after its time");
	}
});

Test test_a4_task1_spline_two_knots_linear_mid("a4.task1.spline.extra.two_knots_mid", []()
{
	Spline<float> s;
	s.set(0.0f, 0.0f);
	s.set(1.0f, 1.0f);
	if (Test::differs(s.at(0.5f), 0.5f))
	{
		throw Test::error("Two-knot mirrored spline should interpolate linearly at midpoint");
	}
});

Test test_a4_task1_cubic_unit_endpoints("a4.task1.spline.extra.cubic_endpoints", []()
{
	Vec3 p0(2.0f, 0.0f, 1.0f);
	Vec3 p1(5.0f, 1.0f, 0.0f);
	Vec3 m0(1.0f, 0.0f, 0.0f);
	Vec3 m1(0.0f, 1.0f, 0.0f);
	Vec3 a0 = Spline<Vec3>::cubic_unit_spline(0.0f, p0, p1, m0, m1);
	Vec3 a1 = Spline<Vec3>::cubic_unit_spline(1.0f, p0, p1, m0, m1);
	if (Test::differs(a0, p0))
	{
		throw Test::error("Hermite spline should match position0 at t=0");
	}
	if (Test::differs(a1, p1))
	{
		throw Test::error("Hermite spline should match position1 at t=1");
	}
});

Test test_a4_task1_spline_collinear_linear("a4.task1.spline.extra.collinear_linear", []()
{
	Spline<float> s;
	s.set(0.0f, 0.0f);
	s.set(1.0f, 1.0f);
	s.set(2.0f, 2.0f);
	for (int i = 1; i <= 9; ++i)
	{
		float t = i * 0.1f;
		float expected = t;
		float actual = s.at(t);
		if (Test::differs(actual, expected))
		{
			throw Test::error(
				"Collinear equally-spaced Catmull-Rom should give linear interpolation at t=" +
				std::to_string(t));
		}
	}
});

Test test_a4_task1_spline_knot_hit_nonuniform("a4.task1.spline.extra.knot_hit_nonuniform", []()
{
	Spline<float> s;
	s.set(0.0f,  1.0f);
	s.set(2.0f,  3.0f);
	s.set(10.0f, 5.0f);
	if (Test::differs(s.at(0.0f),  1.0f)) throw Test::error("Should hit first knot");
	if (Test::differs(s.at(2.0f),  3.0f)) throw Test::error("Should hit middle knot");
	if (Test::differs(s.at(10.0f), 5.0f)) throw Test::error("Should hit last knot");
});
