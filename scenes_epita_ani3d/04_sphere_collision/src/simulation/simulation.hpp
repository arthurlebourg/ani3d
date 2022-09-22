#pragma once

#include "cgp/cgp.hpp"



struct particle_structure
{
    cgp::vec3 p; // Position
    cgp::vec3 v; // Speed

    cgp::vec3 c; // Color
    float r;     // Radius
    float m;     // mass
};

struct plane_structure
{
    cgp::vec3 p;
    cgp::vec3 n;
};

void simulate(std::vector<particle_structure>& particles, std::vector<plane_structure>& walls, float dt);

