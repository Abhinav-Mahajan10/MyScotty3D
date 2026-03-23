
#include "shape.h"
#include "../geometry/util.h"

namespace Shapes {

Vec2 Sphere::uv(Vec3 dir) {
	float u = std::atan2(dir.z, dir.x) / (2.0f * PI_F);
	if (u < 0.0f) u += 1.0f;
	float v = std::acos(-1.0f * std::clamp(dir.y, -1.0f, 1.0f)) / PI_F;
	return Vec2{u, v};
}

BBox Sphere::bbox() const {
	BBox box;
	box.enclose(Vec3(-radius));
	box.enclose(Vec3(radius));
	return box;
}

PT::Trace Sphere::hit(Ray ray) const
{
	//A3T2 - sphere hit
	PT::Trace ret;
	ret.origin = ray.point;
	ret.hit = false;
	ret.distance = 0.0f;
	ret.position = Vec3{};
	ret.normal = Vec3{};
	ret.uv = Vec2{};

	float a = dot(ray.dir, ray.dir);
	float b = 2.0f * dot(ray.point, ray.dir);
	float c = dot(ray.point, ray.point) - radius * radius;
	float disc = b * b - 4.0f * a * c;
	if (disc < 0.0f)
		return ret;

	float root = std::sqrt(disc);
	float inv_2a = 0.5f / a;
	float t0 = (-b - root) * inv_2a;
	float t1 = (-b + root) * inv_2a;
	if (t0 > t1)
		std::swap(t0, t1);

	float t = -1.0f;
	if (t0 >= ray.dist_bounds.x && t0 <= ray.dist_bounds.y)
		t = t0;
	else if (t1 >= ray.dist_bounds.x && t1 <= ray.dist_bounds.y)
		t = t1;
	else
		return ret;

	ret.hit = true;
	ret.distance = t;
	ret.position = ray.at(t);
	ret.normal = ret.position.unit();
	ret.uv = uv(ret.normal);
	return ret;
}

Vec3 Sphere::sample(RNG &rng, Vec3 from) const {
	die("Sampling sphere area lights is not implemented yet.");
}

float Sphere::pdf(Ray ray, Mat4 pdf_T, Mat4 pdf_iT) const {
	die("Sampling sphere area lights is not implemented yet.");
}

Indexed_Mesh Sphere::to_mesh() const {
	return Util::closed_sphere_mesh(radius, 2);
}

} // namespace Shapes

bool operator!=(const Shapes::Sphere& a, const Shapes::Sphere& b) {
	return a.radius != b.radius;
}

bool operator!=(const Shape& a, const Shape& b) {
	if (a.shape.index() != b.shape.index()) return false;
	return std::visit(
		[&](const auto& shape) {
			return shape != std::get<std::decay_t<decltype(shape)>>(b.shape);
		},
		a.shape);
}
