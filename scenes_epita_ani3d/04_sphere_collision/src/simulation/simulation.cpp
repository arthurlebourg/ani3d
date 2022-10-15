#include "simulation.hpp"

using namespace cgp;

Node::Node(cgp::vec3 p, float width)
    : p_(p)
    , center_(p + cgp::vec3(width / 2,width / 2,width / 2))
    , width_(width) 
    , _color_index(int(rand_interval() * 8))
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

bool Node::is_inside_cube(cgp::vec3 b)
{
    if (b.x < p_.x || b.y < p_.y || b.z < p_.z)
    {
        return false;
    }
    if (b.x > p_.x + width_ || b.y > p_.y + width_ || b.z > p_.z + width_)
    {
        return false;
    }
    return true;
}

void Node::add_boule(particle_structure b, std::vector<plane_structure>& walls, float dt, std::shared_ptr<Node> head)
{
    if (!is_inside_cube(b))
    {
        /*std::cout << "ERROR" << std::endl;
        std::cout << "particule: " << b.p << std::endl;
        std::cout << "cube: " << p_ << std::endl;
        std::cout << "width: " << width_ << std::endl << std::endl;
        */
        return;
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
        if (n < boules_per_cube_)
        {   
            boules_.push_back(b);
            simulate(boules_, head, walls, dt);

            if (n == boules_per_cube_ - 1)
            {
                for (auto child : boules_)
                {
                    add_boule(child, walls, dt, head);
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

        children_[i] = std::make_shared<Node>(Node(pos, width_/2));
    }
    children_[i]->add_boule(b, walls, dt, head);
}

std::shared_ptr<Node> Node::get_voisins(cgp::vec3 pos, float width)
{
    if (!is_inside_cube(pos))
    {
        return nullptr;
    }
    //std::cout << "width * 2: " << width * 2 << " wdith_: " << width_ << std::endl;
    size_t i = 0;
    if (pos.x > p_.x + width_/2)
        i+=1;
    if (pos.y > p_.y + width_/2)
        i+=2;
    if (pos.z > p_.z + width_/2)
        i+=4;
    if (width * 2 == width_)
    {
        return children_[i];
    }
    if (children_[i] == nullptr)
    {
        return nullptr;
    }
    return children_[i]->get_voisins(pos, width);
}

std::vector<particle_structure> Node::get_boules(cgp::vec3 pos)
{
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

std::vector<particle_structure> Node::get_boules()
{
    std::vector<particle_structure> res;
    static numarray<vec3> const color_lut = { {0,0,0},{0,0,1},{0,1,0},{0,1,1},{1,0,0},{1,0,1},{1,1,0},{1,1,1} };
    for (auto b : boules_)
    {
        b.c = color_lut[_color_index];
        res.push_back(b);
    }
    for (auto i : children_)
    {
        if (i == nullptr)
        {
            continue;
        }

        auto tmp = i->get_boules();
        for (auto b : tmp)
        {
            res.push_back(b);
        }
    }

    return res;
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

void Node::simulate_opti(float dt, std::vector<plane_structure>& walls, std::vector<particle_structure> &buffer, std::shared_ptr<Node> head)
{
    if (is_leaf())
    {
        simulate(boules_, head, walls, dt);
        std::vector<particle_structure> tmp;
        for (size_t a = 0; a < boules_.size(); a++)
        {
            auto boul = boules_[a];
            if (!is_inside_cube(boul))
            {
                //boules_.erase(boules_.begin() + a);
                buffer.push_back(boul);
            }
            else
            {
                tmp.push_back(boul);
            }
        }
        boules_.clear();
        for (auto i : tmp)
        {
            boules_.push_back(i);
        }
    }
    for (auto i : children_)
    {
        if (i != nullptr)
        {
            i->simulate_opti(dt, walls, buffer, head);
        }
    }
}

void Node::simulate_rec(std::vector<particle_structure>& particles, std::shared_ptr<Node> voisin, std::vector<plane_structure>& walls, float dt)
{
    if (voisin == nullptr)
    {
        return;
    }

    size_t const N = particles.size();
    size_t const M = voisin->boules_.size();
    for (size_t k = 0; k < N; ++k)
	{
		particle_structure& particle = particles[k];
        for (size_t l = 0; l < M; ++l)
        {
            particle_structure& other = voisin->boules_[l];
            collision_boules(particle, other);
        }
    }
    for (auto i : voisin->children_)
    {
        if (i != nullptr)
        {
            simulate_rec(particles, i, walls, dt);
        }
    }
}

void Node::simulate(std::vector<particle_structure>& particles, std::shared_ptr<Node> head, std::vector<plane_structure>& walls, float dt)
{
	size_t const N = particles.size();
    if (N == 0)
    {
        return;
    }

	vec3 const g = { 0,0,-9.81f };

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

    
    for (float i = -width_;  i <= width_; i+= width_)
    {
        for (float j = -width_;  j <= width_; j+= width_)
        {
            for (float h = -width_;  h <= width_; h+= width_)
            {
                if (!(i == 0 && j == 0 && h == 0))
                {
                    auto neighbor = head->get_voisins(vec3(center_.x +i, center_.y + j, center_.z + h), width_);
                    simulate_rec(particles, neighbor, walls, dt);
                }
            }
        }
    }
}
