#pragma once


#include "cgp/cgp.hpp"
#include "environment.hpp"

#include "simulation/simulation.hpp"
#include "implicit_surface/implicit_surface.hpp"

using cgp::mesh_drawable;


struct gui_parameters {
	bool display_frame = true;
	bool add_sphere = true;
	bool move_box = false;
	bool marching_cube = true;
};

// The structure of the custom scene
struct scene_structure : scene_inputs_generic {
	
	// ****************************** //
	// Elements and shapes of the scene
	// ****************************** //
	camera_controller_orbit_euler camera_control;
	camera_projection_perspective camera_projection;
	window_structure window;

	mesh_drawable global_frame;          // The standard global frame
    mesh_drawable implicit_surface; // The shape of the implicit surface
    spatial_domain_grid_3D domain;
	environment_structure environment;   // Standard environment controler
	input_devices inputs;                // Storage for inputs status (mouse, keyboard, window dimension)
	gui_parameters gui;                  // Standard GUI element storage
	cgp::picking_structure picking;
	
	// ****************************** //
	// Elements and shapes of the scene
	// ****************************** //
	cgp::timer_event_periodic timer;
	//std::vector<particle_structure> particles;
	std::shared_ptr<Node> particles;
	std::vector<plane_structure> walls;
	cgp::mesh_drawable sphere;
	cgp::curve_drawable cube_wireframe;

    float sigma = 0.1f;
    float isovalue = 0.05f;


	// ****************************** //
	// Functions
	// ****************************** //

	void initialize();    // Standard initialization to be called before the animation loop
	void display_frame(); // The frame display to be called within the animation loop
	void display_gui();   // The display of the GUI, also called within the animation loop

	void mouse_move_event();
	void mouse_click_event();
	void keyboard_event();
	void idle_frame();

	void emit_particle(std::vector<plane_structure>& walls, float dt);
	void simulation_step(float dt);
	void sphere_display();
};





