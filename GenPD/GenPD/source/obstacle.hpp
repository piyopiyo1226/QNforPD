#pragma once
#include "mesh_struct.hpp"
#include "transformation.hpp"

namespace arcsim {
	struct Obstacle
	{
	public:
		double start_time, end_time;
		bool activated;
		// Gets the last-returned mesh or its transformation
		Mesh &get_mesh();
		const Mesh &get_mesh() const;

		// Gets the state of the mesh at a given time, and updates the internal
		// meshes
		Mesh &get_mesh(double time_sec);

		// lerp with previous mesh at time t - dt
		void blend_with_previous(double t, double dt, double blend);

		const Motion *transform_spline;

		// A mesh containing the original, untransformed object
		Mesh base_mesh;
		// A mesh containing the correct mesh structure
		Mesh curr_state_mesh;

		Obstacle()
			: start_time(0)
			, end_time(0)//end_time(infinity)
			, activated(false) {}


	};
}


