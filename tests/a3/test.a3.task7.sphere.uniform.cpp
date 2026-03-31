#include "test.h"
#include "lib/mathlib.h"
#include "pathtracer/samplers.h"
#include "util/rand.h"

Test test_a3_task7_sphere_uniform_norm("a3.task7.sphere.uniform.norm", []() {
	Samplers::Sphere::Uniform sampler;
	RNG rng(42);
	for (int i = 0; i < 64; i++) {
		Vec3 d = sampler.sample(rng);
		if (!d.valid() || std::abs(d.norm() - 1.0f) > 1e-4f) {
			throw Test::error("Sphere::Uniform::sample should return unit directions!");
		}
		float p = sampler.pdf(d);
		if (!std::isfinite(p) || Test::differs(p, 1.0f / (4.0f * PI_F))) {
			throw Test::error("Sphere::Uniform::pdf should be 1/(4*PI)!");
		}
	}
});
