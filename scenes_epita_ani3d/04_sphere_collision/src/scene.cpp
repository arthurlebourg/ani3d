#include "scene.hpp"


using namespace cgp;

void scene_structure::initialize()
{
	camera_control.initialize(inputs, window); // Give access to the inputs and window global state to the camera controler
	camera_control.set_rotation_axis_z();
	camera_control.look_at({ 3.0f, 2.0f, 2.0f }, {0,0,0}, {0,0,1});
	global_frame.initialize_data_on_gpu(mesh_primitive_frame());


	timer.event_period = 0.5f;

	// Edges of the containing cube in [-1,1]^3
	//  Note: this data structure is set for display purpose - don't use it to compute some information on the cube - it would be un-necessarily complex
	numarray<vec3> border_cube = { {-1,-1,-1},{1,-1,-1}, {1,-1,-1},{1,1,-1}, {1,1,-1},{-1,1,-1}, {-1,1,-1},{-1,-1,-1},
		{-1,-1,1} ,{1,-1,1},  {1,-1,1}, {1,1,1},  {1,1,1}, {-1,1,1},  {-1,1,1}, {-1,-1,1},
		{-1,-1,-1},{-1,-1,1}, {1,-1,-1},{1,-1,1}, {1,1,-1},{1,1,1},   {-1,1,-1},{-1,1,1} };
	cube_wireframe.initialize_data_on_gpu(border_cube);
	cube_wireframe.display_type = curve_drawable_display_type::Segments;

	sphere.initialize_data_on_gpu(mesh_primitive_sphere());

    plane_structure wall1;
    wall1.p = { 0,-1,0 };
    wall1.n = { 0,1,0 };
    walls.push_back(wall1);

    plane_structure wall2;
    wall2.p = { 0,1,0 };
    wall2.n = { 0,-1,0 };
    walls.push_back(wall2);

    plane_structure wall3;
    wall3.p = { -1,0,0 };
    wall3.n = { 1,0,0 };
    walls.push_back(wall3);

    plane_structure wall4;
    wall4.p = { 1,0,0 };
    wall4.n = { -1,0,0 };
    walls.push_back(wall4);

    plane_structure wall5;
    wall5.p = { 0,0,-1 };
    wall5.n = { 0,0,1 };
    walls.push_back(wall5);

    plane_structure wall6;
    wall6.p = { 0,0,1 };
    wall6.n = { 0,0,-1 };
    walls.push_back(wall6);

    int3 const samples = { 25, 25, 25 };
    // Dimension of the domain
    vec3 const length = { 2.5,2.5,2.5 };
    domain = spatial_domain_grid_3D::from_center_length({ 0,0,0 }, length, samples);
	particles = std::make_shared<Node>(Node({-1,-1,-1}, 2));
}

bool first = false;

void scene_structure::display_frame()
{
	// Set the light to the current position of the camera
	environment.light = camera_control.camera_model.position();
	
	if (gui.display_frame)
		draw(global_frame, environment);

	timer.update();

	// Create a new particle if needed
	if (!first)
	{
		emit_particle();
	}	

	// Call the simulation of the particle system
	float const dt = 0.01f * timer.scale;
	//simulate(particles, walls, dt);
	particles->simulate_opti(dt, walls, particles);
	grid_3D<float> field = compute_scalar_field(domain, *particles, sigma);

	// Compute the mesh using marching cube
	mesh m = marching_cube(field, domain, isovalue);
	implicit_surface.initialize_data_on_gpu(m);
    implicit_surface.material.color = {0.83/2,0.94/2,0.97};
    draw(implicit_surface, environment);
	// Display the result
	//sphere_display();
	draw(cube_wireframe, environment);

	if (gui.display_frame)
		draw(global_frame, environment);
}
/*
void scene_structure::sphere_display()
{
	// Display the particles as spheres
	size_t const N = particles.size();
	for (size_t k = 0; k < N; ++k)
	{
		particle_structure const& particle = particles[k];
		sphere.material.color = particle.c;
		sphere.model.translation = particle.p;
		sphere.model.scaling = particle.r;

		draw(sphere, environment);
	}

	// Display the box in which the particles should stay
}
*/
void scene_structure::emit_particle()
{
	// Emit particle with random velocity
	//  Assume first that all particles have the same radius and mass
	static numarray<vec3> const color_lut = { {1,0,0},{0,1,0},{0,0,1},{1,1,0},{1,0,1},{0,1,1} };
	if (timer.event && gui.add_sphere) {
		float const theta = rand_interval(0, 2 * Pi);
		vec3 const v = vec3(1.0f * std::cos(theta), 4.0f, 1.0f * std::sin(theta));

		particle_structure particle;
		particle.p = { 0,0,0 };
		particle.r = 0.08f;//rand_interval(0.08f, 0.16f);
		particle.c = color_lut[int(rand_interval() * color_lut.size())];
		particle.v = v;
		particle.m = rand_interval(1.0f, 2.0f);

		//particles.push_back(particle);
		particles->add_boule(particle);
	}
}


void scene_structure::display_gui()
{
	ImGui::Checkbox("Frame", &gui.display_frame);
	ImGui::SliderFloat("Sigma", &sigma, 0.01f, 2.0f, "%.2f s");
	ImGui::SliderFloat("Isovalue", &isovalue, 0.01f, 1.0f, "%.2f s");
	ImGui::SliderFloat("Time scale", &timer.scale, 0.05f, 2.0f, "%.2f s");
	ImGui::SliderFloat("Time to add new sphere", &timer.event_period, 0.05f, 2.0f, "%.2f s");
	ImGui::Checkbox("Add sphere", &gui.add_sphere);
	ImGui::Checkbox("Move box", &gui.move_box);
}

void scene_structure::mouse_move_event()
{
    if (gui.move_box) 
    {
        vec2 const& p = inputs.mouse.position.current;

		if (inputs.mouse.click.left)
        {
			vec2 const translation_screen = p - picking.screen_clicked;
            float val = translation_screen.y * 0.1;
            cube_wireframe.model.rotation *= rotation_transform::from_axis_angle({1,0,0}, val);
            for (size_t i = 0; i < 6; i++)
            {
                walls[i].p = rotation_around_center(rotation_transform::from_axis_angle({1,0,0}, val), {0,0,0}) * walls[i].p;

                walls[i].n = rotation_transform::from_axis_angle({1,0,0}, val) * walls[i].n;

            }

        }
    }
    else
    {
        camera_control.action_mouse_move(environment.camera_view);
    }
}
void scene_structure::mouse_click_event()
{
	camera_control.action_mouse_click(environment.camera_view);
	if (inputs.mouse.click.last_action == last_mouse_cursor_action::release_left)
	{
		// Releasing the left click indicates the end of the deformation: disable the picking, save the new position and update the normals
		picking.active = false;
	}
}
void scene_structure::keyboard_event()
{
	camera_control.action_keyboard(environment.camera_view);
}
void scene_structure::idle_frame()
{
	camera_control.idle_frame(environment.camera_view);
}

