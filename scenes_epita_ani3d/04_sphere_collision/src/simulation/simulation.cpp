#include "simulation.hpp"

using namespace cgp;

void collision_boules(particle_structure &boule1, particle_structure &boule2)
{
    float norm_boules = norm(boule2.p - boule1.p);
    if (norm_boules <= boule1.r + boule2.r)
    {
        cgp::vec3 u = (boule1.p - boule2.p) / norm_boules;

        if (norm(boule1.v - boule2.v) > 0.5)
        {
            float j = 2 * dot(boule2.v - boule1.v, u) / (1/boule1.m + 1/boule2.m); 
            boule1.v = boule1.v + j/boule1.m * u;
            boule2.v = boule2.v - j/boule2.m * u;
        }
        else
        {
            boule1.v = 0.5*boule1.v;
            boule2.v = 0.5*boule2.v;
        }

        float d = boule1.r + boule2.r  - norm_boules;
        boule1.p += d/2*u;
        boule2.p -= d/2*u;
    }
}


void simulate(std::vector<particle_structure>& particles, std::vector<plane_structure>& walls, float dt)
{
	vec3 const g = { 0,0,-9.81f };
	size_t const N = particles.size();
	for (size_t k = 0; k < N; ++k)
	{
		particle_structure& particle = particles[k];

		vec3 const f = particle.m * g;

		particle.v = (1 - 0.9f * dt) * particle.v + dt * f;
		particle.p = particle.p + dt * particle.v;
	}

	// To do :
	//  Handle collision ...
	for (size_t k = 0; k < N; ++k)
	{
		particle_structure& particle = particles[k];
        for (size_t i = 0; i < 6; ++i)
        {
            plane_structure& wall = walls[i];
            float detection = dot(particle.p - wall.p, wall.n);
            if (detection <= particle.r)
            {
                float friction = 1;
                float impact = 1;
                auto tmp = dot(particle.v, wall.n) * wall.n;
                vec3 new_v = friction * (particle.v - tmp) - impact * tmp;
                particle.v = 0.9f*new_v;

                float d = particle.r - dot(particle.p - wall.p, wall.n);
                particle.p += d * wall.n;
            }
        }

        for (size_t i = 0; i < k; ++i)
        {
            particle_structure& other = particles[i];
            collision_boules(particle, other);
        }
        for (size_t i = k+1; i < N; ++i)
        {
            particle_structure& other = particles[i];
            collision_boules(particle, other);
        }

    }

}
