#include "arcsim_simulation.hpp"
#include "arcsim_handle.hpp"
#include "remesh.hpp"
#include "collision.hpp"
#include "near_obs.hpp"

namespace arcsim {
	void remeshing_step(Arc_sim &sim)//���Ƃ�sim�@�̂��ɂ܂킵�Ă���
	{//
	 // remesh
	 //----------------------------------------------------------------------------//
		for (size_t c = 0; c < sim.cloths.size(); c++) {

			std::vector<AccelStruct *> obs_accs = create_accel_structs(sim.obstacle_meshes, false);//accel thestructure
			std::map<Node *, Plane> planes = nearest_obstacle_planes(sim.cloths[c].mesh.nodes, obs_accs);																					  
			destroy_accel_structs(obs_accs);
			dynamic_remesh(sim.cloths[c].mesh, planes);
		}
		//--------------------------------------------------------collision--------//
	}
	//---------------------------------------in upper------------------------
	//dynamic_remesh(sim.cloths[c].mesh, planes);
	//------------------------------------------------------------------------

	// separate
	//	if (sim.enabled[separation]) {
	//		separate(sim.cloth_meshes, sim.obstacle_meshes);
	//	}
	// apply pop filter ������g���Ă��Ȃ���x�n�@������ŏ��͂Ƃ΂����Ƃɂ���ۂ񂿂܂�

	//if (sim.enabled[popfilter] && !initializing) {
	//	vector<Constraint *> cons = get_constraints(sim, true);
	//	for (size_t c = 0; c < sim.cloths.size(); c++) apply_pop_filter(sim.cloths[c], cons);
	//	delete_constraints(cons);
	//}

// pop fileter �K�i�j}
}