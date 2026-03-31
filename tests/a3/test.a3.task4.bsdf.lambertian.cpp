#include "test.h"
#include "scene/material.h"
#include "scene/texture.h"
#include "util/rand.h"

Test test_a3_task4_bsdf_lambertian_simple("a3.task4.bsdf.lambertian.simple", []() {
	// This test just checks that the sample function produces a valid sample.

	auto alb_t = std::make_shared<Texture>(Texture{Textures::Constant{Spectrum{1.0f}}});
	auto bsdf = Materials::Lambertian{alb_t};

	Vec3 out;

	RNG rng(1);

	Materials::Scatter s = bsdf.scatter(rng, out, {});
	// Check that the direction is valid
	if (!s.direction.valid() || s.direction.norm() == 0.0f) {
		throw Test::error("BSDF produced invalid sample!");
	}
	float pdf = bsdf.pdf(out, s.direction);
	// Check that the pdf is valid
	if (!std::isfinite(pdf) || pdf < 0.0f) {
		throw Test::error("BSDF produced sample with invalid pdf!");
	}
	// Check the value against the exact attenuation for the first sample from RNG with seed 1
	if (Test::differs(s.attenuation, Spectrum{0.317861f, 0.317861f, 0.317861f})) {
		throw Test::error("BSDF sample attenuation was not equivalent to evaluate!");
	}
});

Test test_a3_task4_bsdf_lambertian_pdf_lower_hemisphere("a3.task4.bsdf.lambertian.pdf.lower_hemisphere", []() {
	auto alb_t = std::make_shared<Texture>(Texture{Textures::Constant{Spectrum{1.0f}}});
	auto bsdf = Materials::Lambertian{alb_t};
	Vec3 out{};
	Vec3 in_down = Vec3{0.1f, -1.0f, 0.1f}.unit();
	float pdf = bsdf.pdf(out, in_down);
	if (pdf != 0.0f) {
		throw Test::error("Cosine hemisphere PDF should be zero for directions below the surface!");
	}
	Spectrum ev = bsdf.evaluate(out, in_down, {});
	if (ev.luma() != 0.0f) {
		throw Test::error("Lambertian evaluate should be zero when cos(theta) <= 0!");
	}
});

Test test_a3_task4_bsdf_lambertian_evaluate_white("a3.task4.bsdf.lambertian.evaluate.grazing", []() {
	auto alb_t = std::make_shared<Texture>(Texture{Textures::Constant{Spectrum{1.0f}}});
	auto bsdf = Materials::Lambertian{alb_t};
	Vec3 out{};
	Vec3 in_up = Vec3{0.f, 1.f, 0.f};
	Spectrum ev = bsdf.evaluate(out, in_up, {});
	if (Test::differs(ev, Spectrum{1.0f / PI_F, 1.0f / PI_F, 1.0f / PI_F})) {
		throw Test::error("Lambertian at normal incidence with white albedo should be albedo/PI!");
	}
});
