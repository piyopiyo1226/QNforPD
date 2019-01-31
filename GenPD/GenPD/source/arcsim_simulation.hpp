#pragma once
#include "obstacle.hpp"
#include "arcsim_handle.hpp"



namespace arcsim {
	struct Arc_sim
	{
		double time;
		int frame, step;

		std::vector<Cloth> cloths;

		int frame_step, save_every;
		double frame_time, step_time;
		double end_time, end_frame;
		double passive_time;

		std::vector<Handle * > handles;
		std::vector<Obstacle>  obstacles;
		//	std::vector<Morph> morphs;

		Vec3 gravity;
		//	Wind wind;
		double friction, obs_friction;

		enum {
			Proximity,
			Physics,
			StarinLimiting,
			Collision,
			Remeshing,
			Separation,
			PopFilter,
			Plasticity,
			Fracture,
			nModelues
		};

		bool enabled[nModelues];

		//		Timer timers[nModules];
		std::vector<Mesh * > cloth_meshes, obstacle_meshes;
	};

	//namespace global_index {
		extern Arc_sim sim;
	//}
}