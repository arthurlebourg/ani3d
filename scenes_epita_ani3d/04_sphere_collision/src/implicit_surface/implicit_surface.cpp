#include "implicit_surface.hpp"


using namespace cgp;

float field_evaluate(vec3 const& p, std::vector<particle_structure> particles);

grid_3D<float> compute_scalar_field(spatial_domain_grid_3D const& domain,std::vector<particle_structure> particles)
{
	grid_3D<float> field;
	field.resize(domain.samples);

	// Fill the discrete field values
	for (int kz = 0; kz < domain.samples.z; kz++) {
		for (int ky = 0; ky < domain.samples.y; ky++) {
			for (int kx = 0; kx < domain.samples.x; kx++) {

				vec3 const p = domain.position({ kx, ky, kz });
				field(kx, ky, kz) = field_evaluate(p, particles);

			}
		}
	}

	return field;
}


// Parameterization of a basic "blob"-function - Gaussian centered at point p0
float blob(vec3 const& p, vec3 const& p0)
{
    float sigma = 1;
	float const d = norm(p - p0);
	float const value = std::exp(-(d * d)/(sigma*sigma));
	return value;
}

// Evaluate the field at an arbitrary position in space
float field_evaluate(vec3 const& p, std::vector<particle_structure> particles)
{
	float value = 0.0f;
	size_t const N = particles.size();
	for (size_t k = 0; k < N; ++k)
	{
		particle_structure const& particle = particles[k];
        value = std::max(value, blob(p,particle.p));

	}

	return value;
}
