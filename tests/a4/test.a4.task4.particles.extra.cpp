#include "test.h"
#include "geometry/util.h"
#include "lib/mathlib.h"
#include "pathtracer/aggregate.h"
#include "scene/particles.h"

Test test_a4_task4_particles_lifetime_expires("a4.task4.particles.lifetime.expires", []() {
	PT::Aggregate empty;
	Particles particles;
	particles.particles = {{Vec3(0.0f), Vec3(0.0f), 0.0005f}};
	bool alive = particles.particles[0].update(empty, particles.gravity, particles.radius, 0.01f);
	if (alive)
	{
		throw Test::error("Particle should be removed when age falls through zero.");
	}
});

Test test_a4_task4_particles_resting_slides_down("a4.task4.particles.resting_accelerates", []() {
	PT::Aggregate empty;
	Particles particles;
	particles.particles = {{Vec3(0.0f, 2.0f, 0.0f), Vec3(0.0f), 10.0f}};
	Particles::Particle p = particles.particles[0];
	bool alive = p.update(empty, particles.gravity, particles.radius, 0.01f);
	if (!alive)
	{
		throw Test::error("Unexpected particle death.");
	}
	if (p.velocity.y >= -1e-6f)
	{
		throw Test::error("Stationary particle should accelerate downward under gravity.");
	}
});

Test test_a4_task4_particles_reflection_upward(
	"a4.task4.particles.reflection.upward", []()
	{
		PT::Tri_Mesh ground_mesh{Util::square_mesh(10), false};
		PT::Aggregate ground{PT::List{std::vector{PT::Instance{&ground_mesh, nullptr, Mat4::I}}}};

		Particles::Particle p;
		p.position = Vec3(0.0f, 1.0f, 0.0f);
		p.velocity = Vec3(0.0f, -3.0f, 0.0f);
		p.age = 100.0f;

		Vec3 gravity(0.0f, -9.8f, 0.0f);
		float radius = 0.2f;

		bool bounced = false;
		for (int i = 0; i < 500 && p.update(ground, gravity, radius, 0.01f); ++i)
		{
			if (p.velocity.y > 1e-4f)
			{
				bounced = true;
				break;
			}
		}
		if (!bounced)
		{
			throw Test::error("Particle should bounce off the floor with upward velocity.");
		}
	});

Test test_a4_task4_particles_no_gravity_straight_line(
	"a4.task4.particles.no_gravity.straight_line", []()
	{
		PT::Aggregate empty;
		Vec3 gravity(0.0f, 0.0f, 0.0f);
		float radius = 0.01f;

		Particles::Particle p;
		p.position = Vec3(0.0f, 0.0f, 0.0f);
		p.velocity = Vec3(1.0f, 2.0f, 3.0f);
		p.age = 100.0f;

		float dt = 0.05f;
		for (int i = 0; i < 10; ++i)
		{
			Vec3 expected = p.position + p.velocity * dt;
			p.update(empty, gravity, radius, dt);
			if (Test::differs(p.position, expected))
			{
				throw Test::error("Zero-gravity free particle must move in a straight line.");
			}
		}
	});
