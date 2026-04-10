
#include "../geometry/spline.h"

template<typename T> T Spline<T>::at(float time) const
{

	// A4T1b: Evaluate a Catumull-Rom spline

	// Given a time, find the nearest positions & tangent values
	// defined by the control point map.

	// Transform them for use with cubic_unit_spline

	// Be wary of edge cases! What if time is before the first knot,
	// before the second knot, etc...

	if (knots.empty())
	{
		return T();
	}

	auto first = knots.begin();
	if (knots.size() == 1)
	{
		return first->second;
	}

	auto last = knots.rbegin();
	if (time <= first->first)
	{
		return first->second;
	}
	if (time >= last->first)
	{
		return last->second;
	}

	auto k2_it = knots.upper_bound(time);
	auto k1_it = std::prev(k2_it);
	float t1 = k1_it->first;
	float t2 = k2_it->first;
	const T &p1 = k1_it->second;
	const T &p2 = k2_it->second;

	T p0;
	float t0;
	if (k1_it == first)
	{
		t0 = t1 - (t2 - t1);
		p0 = p1 - (p2 - p1);
	}
	else
	{
		auto k0_it = std::prev(k1_it);
		t0 = k0_it->first;
		p0 = k0_it->second;
	}

	T p3;
	float t3;
	if (k2_it == std::prev(knots.end()))
	{
		t3 = t2 + (t2 - t1);
		p3 = p2 + (p2 - p1);
	}
	else
	{
		auto k3_it = std::next(k2_it);
		t3 = k3_it->first;
		p3 = k3_it->second;
	}

	T m0 = (p2 - p0) / (t2 - t0);
	T m1 = (p3 - p1) / (t3 - t1);
	float span = t2 - t1;
	float u = (time - t1) / span;

	return cubic_unit_spline(u, p1, p2, m0 * span, m1 * span);
}

template<typename T>
T Spline<T>::cubic_unit_spline(float time, const T& position0, const T& position1,
                               const T& tangent0, const T& tangent1)
{

	// A4T1a: Hermite Curve over the unit interval

	// Given time in [0,1] compute the cubic spline coefficients and use them to compute
	// the interpolated value at time 'time' based on the positions & tangents

	// Note that Spline is parameterized on type T, which allows us to create splines over
	// any type that supports the * and + operators.

	float t = time;
	float t2 = t * t;
	float t3 = t2 * t;

	float h00 = 2.0f * t3 - 3.0f * t2 + 1.0f;
	float h10 = t3 - 2.0f * t2 + t;
	float h01 = -2.0f * t3 + 3.0f * t2;
	float h11 = t3 - t2;

	return h00 * position0 + h10 * tangent0 + h01 * position1 + h11 * tangent1;
}

template class Spline<float>;
template class Spline<double>;
template class Spline<Vec4>;
template class Spline<Vec3>;
template class Spline<Vec2>;
template class Spline<Mat4>;
template class Spline<Spectrum>;
