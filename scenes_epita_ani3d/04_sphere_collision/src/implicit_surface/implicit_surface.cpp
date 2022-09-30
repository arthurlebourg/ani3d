#include "implicit_surface.hpp"


using namespace cgp;

float field_evaluate(vec3 const& p, Node particles, const float &sigma);

grid_3D<float> compute_scalar_field(spatial_domain_grid_3D const& domain, Node particles,
const float &sigma)
{
	grid_3D<float> field;
	field.resize(domain.samples);

	// Fill the discrete field values
	for (int kz = 0; kz < domain.samples.z; kz++) {
		for (int ky = 0; ky < domain.samples.y; ky++) {
			for (int kx = 0; kx < domain.samples.x; kx++) {

				vec3 const p = domain.position({ kx, ky, kz });
				field(kx, ky, kz) = field_evaluate(p, particles, sigma);

			}
		}
	}

	return field;
}


// Parameterization of a basic "blob"-function - Gaussian centered at point p0
float blob(vec3 const& p, vec3 const& p0, const float &sigma)
{
	float const d = norm(p - p0);
	//float const value = 1 - sqrt(d);
	float const value = exp(-(d*d)/(sigma*sigma));
	return value;
}

// Evaluate the field at an arbitrary position in space
float field_evaluate(vec3 const& p, Node particles, const float &sigma)
{
	float value = 0.0f;
	auto particles_vec = particles.get_boules(p);
	size_t const N = particles_vec.size();
	for (size_t k = 0; k < N; ++k)
	{
		particle_structure const& particle = particles_vec[k];
        value = std::max(value, blob(p,particle.p, sigma));
	}
	return value;
}
