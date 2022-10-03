#include "simulation.hpp"

using namespace cgp;

Node::Node(cgp::vec3 p, float width)
    : p_(p)
    , width_(width) 
{
    children_ = std::vector<std::shared_ptr<Node>>(8, nullptr);
}

bool Node::is_leaf()
{
    for (auto i : children_)
    {
        if (i != nullptr)
        {
            return false;
        }
    }
    return true;
}

std::vector<std::shared_ptr<Node>> Node::get_leaves()
{
    std::vector<std::shared_ptr<Node>> nodes;
    for (auto i : children_)
    {
        if (i == nullptr)
            continue;
        if (i->is_leaf())
        {
            nodes.push_back(i);
        }
        else
        {
            for (auto j : i->get_leaves())
                nodes.push_back(j);
        }
    }

    return nodes;
}

bool Node::is_inside_cube(particle_structure b)
{
    if (b.p.x < p_.x || b.p.y < p_.y || b.p.z < p_.z)
    {
        return false;
    }
    if (b.p.x > p_.x + width_ || b.p.y > p_.y + width_ || b.p.z > p_.z + width_)
    {
        return false;
    }
    return true;
}

void Node::add_boule(particle_structure b)
{
    if (!is_inside_cube(b))
    {
        /*std::cout << "ERROR" << std::endl;
        std::cout << "particule: " << b.p << std::endl;
        std::cout << "cube: " << p_ << std::endl;
        std::cout << "width: " << width_ << std::endl;*/
    }
    size_t i = 0;
    if (b.p.x > p_.x + width_/2)
        i+=1;
    if (b.p.y > p_.y + width_/2)
        i+=2;
    if (b.p.z > p_.z + width_/2)
        i+=4;

    if (is_leaf())
    {
        size_t n = boules_.size();
        if (width_ >= 0.5)
        {
            std::cout << "cube: " << p_ << std::endl;
            std::cout << "boules_ size: " << n << std::endl;
            std::cout << "width: " << width_ << std::endl;
            std::cout << "particle: " << b.p << std::endl << std::endl;
        }
        if (n < 8)
        {   
            boules_.push_back(b);
            if (n == 7)
            {
                for (auto child : boules_)
                {
                    add_boule(child);
                }
                boules_.clear();
            }
            return;
        }
    }
    if (children_[i] == nullptr)
    {
        cgp::vec3 pos = p_;
        if (b.p.x > p_.x + width_/2)
            pos.x += width_/2;
        if (b.p.y > p_.y + width_/2)
            pos.y += width_/2;
        if (b.p.z > p_.z + width_/2)
            pos.z += width_/2;
        if (width_ != 0 && false)
        {
            std::cout << "cube: " << p_ << std::endl;
            std::cout << "width: " << width_ << std::endl;
            std::cout << "particle: " << b.p << std::endl;
            std::cout << "pos: " << pos << std::endl << std::endl;
        }
        children_[i] = std::make_shared<Node>(Node(pos, width_/2));
    }
    children_[i]->add_boule(b);
}

std::vector<particle_structure> Node::get_boules(cgp::vec3 pos)
{
    if (is_leaf())
    {
        return boules_;
    }
    size_t i = 0;
    if (pos.x > p_.x + width_/2)
        i+=1;
    if (pos.y > p_.y + width_/2)
        i+=2;
    if (pos.z > p_.z + width_/2)
        i+=4;
    
    if (children_[i] == nullptr)
    {
        return boules_;
    }

    return children_[i]->get_boules(pos);
}

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

void Node::simulate_opti(float dt, std::vector<plane_structure>& walls, std::shared_ptr<Node> head)
{
    for (auto i : children_)
    {
        if (i == nullptr)
        {
            continue;
        }

        if (i->is_leaf())
        {
            simulate(i->boules_, walls, dt);
            for (auto child : boules_)
            {
                head->add_boule(child);
            }
            i->children_.empty();
            // move boules in the octree, their location changed
        }
        else
        {
            i->simulate_opti(dt, walls, head);
        }
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
                float friction = 0.2;
                float impact = 0.5;
                auto tmp = dot(particle.v, wall.n) * wall.n;
                vec3 new_v = friction * (particle.v - tmp) - impact * tmp;
                particle.v = 0.5f*new_v;

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
