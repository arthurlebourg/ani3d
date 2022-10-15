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


class Node
{
public:
    Node(cgp::vec3 p, float width);

    bool is_leaf();

    std::vector<std::shared_ptr<Node>> get_leaves();

    bool is_inside_cube(particle_structure b);
    bool is_inside_cube(cgp::vec3 b);


    void add_boule(particle_structure b, std::vector<plane_structure>& walls, float dt, std::shared_ptr<Node> head);

    std::shared_ptr<Node> get_voisins(cgp::vec3 pos, float width);

    std::vector<particle_structure> get_boules(cgp::vec3 pos);
    std::vector<particle_structure> get_boules();


    void simulate_opti(float dt, std::vector<plane_structure>& walls, std::vector<particle_structure> &buffer, std::shared_ptr<Node> head);
    void simulate_rec(std::vector<particle_structure>& particles, std::shared_ptr<Node> voisin, std::vector<plane_structure>& walls, float dt);
    void simulate(std::vector<particle_structure>& particles, std::shared_ptr<Node> head, std::vector<plane_structure>& walls, float dt);


private:
    cgp::vec3 p_;
    cgp::vec3 center_;
    float width_;
    std::vector<particle_structure> boules_;
    std::vector<std::shared_ptr<Node>> children_;
    cgp::vec3 color_;
    

    size_t boules_per_cube_ = 20;
};

void simulate(std::vector<particle_structure>& particles, std::shared_ptr<Node> head, float width, std::vector<plane_structure>& walls, float dt);




