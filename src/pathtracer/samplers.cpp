
#include "samplers.h"
#include "../scene/shape.h"
#include "../util/rand.h"

#include <algorithm>
#include <cmath>
#include <deque>
#include <vector>

constexpr bool IMPORTANCE_SAMPLING = true;
constexpr bool USE_ALIAS_ENV_SAMPLER = true;

namespace {

static void build_vose_alias_table(const std::vector<float>& weights, std::vector<float>& prob,
                                   std::vector<uint32_t>& alt)
{
	const size_t n = weights.size();
	prob.assign(n, 0.f);
	alt.assign(n, 0u);
	if (n == 0)
	{
		return;
	}
	float sum = 0.f;
	for (float w : weights) sum += w;
	if (sum <= 1e-30f)
	{
		prob.assign(n, 1.f / float(n));
		for (size_t i = 0; i < n; i++) alt[i] = 0u;
		return;
	}
	std::vector<float> scaled(n);
	float scale = float(n) / sum;
	for (size_t i = 0; i < n; i++) scaled[i] = weights[i] * scale;
	std::deque<size_t> small, large;
	for (size_t i = 0; i < n; i++)
	{
		if (scaled[i] < 1.f)
		{
			small.push_back(i);
		}
		else
		{
			large.push_back(i);
		}
	}
	while (!small.empty() && !large.empty())
	{
		size_t s = small.front();
		small.pop_front();
		size_t l = large.front();
		large.pop_front();
		prob[s] = scaled[s];
		alt[s] = static_cast<uint32_t>(l);
		scaled[l] = scaled[l] + scaled[s] - 1.f;
		if (scaled[l] < 1.f)
		{
			small.push_back(l);
		}
		else
		{
			large.push_back(l);
		}
	}
	while (!large.empty())
	{
		size_t l = large.front();
		large.pop_front();
		prob[l] = 1.f;
	}
	while (!small.empty())
	{
		size_t s = small.front();
		small.pop_front();
		prob[s] = 1.f;
	}
}

} // namespace

namespace Samplers {

Vec2 Rect::sample(RNG &rng) const
{
	//A3T1 - step 2 - supersampling
	return Vec2(rng.unit() * size.x, rng.unit() * size.y);
}

float Rect::pdf(Vec2 at) const {
	if (at.x < 0.0f || at.x > size.x || at.y < 0.0f || at.y > size.y) return 0.0f;
	return 1.0f / (size.x * size.y);
}

Vec2 Circle::sample(RNG &rng) const
{
	//A3EC - bokeh - circle sampling
	float r = radius * std::sqrt(rng.unit());
	float theta = 2.0f * PI_F * rng.unit();
	return center + Vec2(r * std::cos(theta), r * std::sin(theta));
}

float Circle::pdf(Vec2 at) const
{
	//A3EC - bokeh - circle pdf
	Vec2 d = at - center;
	if (dot(d, d) > radius * radius)
		return 0.0f;
	return 1.0f / (PI_F * radius * radius);
}

Vec3 Point::sample(RNG &rng) const {
	return point;
}

float Point::pdf(Vec3 at) const {
	return at == point ? 1.0f : 0.0f;
}

Vec3 Triangle::sample(RNG &rng) const {
	float u = std::sqrt(rng.unit());
	float v = rng.unit();
	float a = u * (1.0f - v);
	float b = u * v;
	return a * v0 + b * v1 + (1.0f - a - b) * v2;
}

float Triangle::pdf(Vec3 at) const {
	float a = 0.5f * cross(v1 - v0, v2 - v0).norm();
	float u = 0.5f * cross(at - v1, at - v2).norm() / a;
	float v = 0.5f * cross(at - v2, at - v0).norm() / a;
	float w = 1.0f - u - v;
	if (u < 0.0f || v < 0.0f || w < 0.0f) return 0.0f;
	if (u > 1.0f || v > 1.0f || w > 1.0f) return 0.0f;
	return 1.0f / a;
}

Vec3 Hemisphere::Uniform::sample(RNG &rng) const {

	float Xi1 = rng.unit();
	float Xi2 = rng.unit();

	float theta = std::acos(Xi1);
	float phi = 2.0f * PI_F * Xi2;

	float xs = std::sin(theta) * std::cos(phi);
	float ys = std::cos(theta);
	float zs = std::sin(theta) * std::sin(phi);

	return Vec3(xs, ys, zs);
}

float Hemisphere::Uniform::pdf(Vec3 dir) const {
	if (dir.y < 0.0f) return 0.0f;
	return 1.0f / (2.0f * PI_F);
}

Vec3 Hemisphere::Cosine::sample(RNG &rng) const {

	float phi = rng.unit() * 2.0f * PI_F;
	float cos_t = std::sqrt(rng.unit());

	float sin_t = std::sqrt(1 - cos_t * cos_t);
	float x = std::cos(phi) * sin_t;
	float z = std::sin(phi) * sin_t;
	float y = cos_t;

	return Vec3(x, y, z);
}

float Hemisphere::Cosine::pdf(Vec3 dir) const {
	if (dir.y < 0.0f) return 0.0f;
	return dir.y / PI_F;
}

Vec3 Sphere::Uniform::sample(RNG &rng) const
{
	//A3T7 - sphere sampler

	// Generate a uniformly random point on the unit sphere.
	// Tip: start with Hemisphere::Uniform

	float z = 1.f - 2.f * rng.unit();
	float phi = 2.f * PI_F * rng.unit();
	float r = std::sqrt(std::max(0.f, 1.f - z * z));
	return Vec3(r * std::cos(phi), z, r * std::sin(phi)).unit();
}

float Sphere::Uniform::pdf(Vec3 dir) const
{
	return 1.0f / (4.0f * PI_F);
}

Sphere::Image::Image(const HDR_Image& image)
{
	//A3T7 - image sampler init

	// Set up importance sampling data structures for a spherical environment map image.
	// You may make use of the _pdf, _cdf, and total members, or create your own.

	const auto [_w, _h] = image.dimension();
	w = _w;
	h = _h;
	const size_t n = size_t(w) * size_t(h);
	if (n == 0)
	{
		return;
	}
	std::vector<float> weights(n, 0.f);
	for (uint32_t y = 0; y < h; y++)
	{
		float v_tex = (float(y) + 0.5f) / float(h);
		float sin_t = std::sin(PI_F * v_tex);
		for (uint32_t x = 0; x < w; x++)
		{
			float lu = image.at(x, y).luma();
			weights[size_t(y) * w + x] = lu * sin_t;
		}
	}
	float sum = 0.f;
	for (float t : weights) sum += t;
	if (sum <= 1e-30f)
	{
		float u = 1.f / float(n);
		_pdf.assign(n, u);
	}
	else
	{
		_pdf.resize(n);
		for (size_t i = 0; i < n; i++)
		{
			_pdf[i] = weights[i] / sum;
		}
	}
	_cdf.resize(n);
	float c = 0.f;
	for (size_t i = 0; i < n; i++)
	{
		c += _pdf[i];
		_cdf[i] = c;
	}
	if (!_cdf.empty())
	{
		_cdf.back() = 1.f;
	}
	build_vose_alias_table(weights, _alias_prob, _alias_alt);
}

static Vec3 direction_from_texel(uint32_t xi, uint32_t yi, uint32_t w, uint32_t h)
{
	float u_tex = (float(xi) + 0.5f) / float(w);
	float v_tex = (float(yi) + 0.5f) / float(h);
	float phi = u_tex * 2.f * PI_F;
	float dir_y = -std::cos(PI_F * v_tex);
	float sin_t = std::sin(PI_F * v_tex);
	return Vec3(std::cos(phi) * sin_t, dir_y, std::sin(phi) * sin_t).unit();
}

Vec3 Sphere::Image::sample(RNG &rng) const
{
	Sphere::Uniform uni;
	if (!IMPORTANCE_SAMPLING)
	{
		// Step 1: Uniform sampling
		// Declare a uniform sampler and return its sample
		return uni.sample(rng);
	}
	const size_t n = size_t(w) * size_t(h);
	if (n == 0)
	{
		return uni.sample(rng);
	}
	size_t idx = 0;
	// Step 2: Importance sampling
	// Use your importance sampling data structure to generate a sample direction.
	// Tip: std::upper_bound
	if constexpr (USE_ALIAS_ENV_SAMPLER)
	{
		if (!_alias_prob.empty())
		{
			uint32_t i = rng.integer(0, int32_t(n));
			float u = rng.unit();
			idx = (u < _alias_prob[i]) ? size_t(i) : size_t(_alias_alt[i]);
		}
		else
		{
			float u = rng.unit();
			auto it = std::upper_bound(_cdf.begin(), _cdf.end(), u);
			idx = std::min<size_t>(size_t(it - _cdf.begin()), n - 1);
		}
	}
	else
	{
		float u = rng.unit();
		auto it = std::upper_bound(_cdf.begin(), _cdf.end(), u);
		idx = std::min<size_t>(size_t(it - _cdf.begin()), n - 1);
	}
	uint32_t yi = uint32_t(idx / w);
	uint32_t xi = uint32_t(idx % w);
	return direction_from_texel(xi, yi, w, h);
}

float Sphere::Image::pdf(Vec3 dir) const
{
	Sphere::Uniform uni;
	if (!IMPORTANCE_SAMPLING)
	{
		// Step 1: Uniform sampling
		// Declare a uniform sampler and return its pdf
		return uni.pdf(dir);
	}
	if (w == 0 || h == 0 || _pdf.empty())
	{
		return uni.pdf(dir);
	}
	// A3T7 - image sampler importance sampling pdf
	// What is the PDF of this distribution at a particular direction?
	Vec2 uv = Shapes::Sphere::uv(dir);
	float u = uv.x;
	float v = uv.y;
	if (u >= 1.f)
	{
		u = std::nextafter(1.f, 0.f);
	}
	if (v >= 1.f)
	{
		v = std::nextafter(1.f, 0.f);
	}
	uint32_t xi = std::min(w - 1u, uint32_t(u * float(w)));
	uint32_t yi = std::min(h - 1u, uint32_t(v * float(h)));
	size_t cell = size_t(yi) * w + xi;
	float p_pixel = _pdf[cell];
	float sin_theta = std::sqrt(std::max(0.f, 1.f - dir.y * dir.y));
	sin_theta = std::max(sin_theta, 1e-8f);
	float jacobian = (float(w) * float(h)) / (2.f * PI_F * PI_F * sin_theta);
	return p_pixel * jacobian;
}

} // namespace Samplers
