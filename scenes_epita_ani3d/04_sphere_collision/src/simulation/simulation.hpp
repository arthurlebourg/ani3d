#pragma once

#include "cgp/cgp.hpp"
#include <memory>
#define MAX_DEPTH 8

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

class Node
{
public:
    Node(cgp::vec3 p, float width);

    bool is_leaf();

    std::vector<std::shared_ptr<Node>> get_leaves();

    bool is_inside_cube(particle_structure b);

    void add_boule(particle_structure b);

    std::vector<particle_structure> get_boules(cgp::vec3 pos);
    std::vector<particle_structure> get_boules();


    void simulate_opti(float dt, std::vector<plane_structure>& walls, std::shared_ptr<Node> head);

private:
    cgp::vec3 p_;
    float width_;
    std::vector<particle_structure> boules_;
    std::vector<std::shared_ptr<Node>> children_;

    size_t boules_per_cube_ = 20;
};





