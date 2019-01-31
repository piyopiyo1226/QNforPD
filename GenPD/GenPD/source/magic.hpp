#ifndef MAGIC_HPP
#define MAGIC_HPP
namespace arcsim {

	struct Magic {//Ç±Ç±ÇÕÇ÷ÇÒÇ∑Ç§Ç™ÇΩÇ≠Ç≥ÇÒÇÕÇ¢Ç¡ÇƒÇ§ÇíÅ@
		bool fixed_high_res_mesh;
		double handle_stiffness, collision_stiffness;
		double repulsion_thickness, projection_thickness;
		double edge_flip_threshold;
		double rib_stiffening;
		bool combine_tensors;
		bool preserve_creases;
		bool add_jitter;
		double separation_step_size;
		int relax_method, max_cracks;
		bool enable_localopt;
		Magic()
			: fixed_high_res_mesh(false)
			, handle_stiffness(1e3)
			, collision_stiffness(1e9)
			, repulsion_thickness(1e-3)
			, projection_thickness(1e-4)
			, edge_flip_threshold(1e-2)
			, rib_stiffening(1)
			, combine_tensors(true)
			, preserve_creases(false)
			, add_jitter(false)
			, separation_step_size(1e-2)
			, relax_method(0)
			, max_cracks(100)
			, enable_localopt(false) {}
	};

	namespace global_index {
		extern Magic magic;
	}
}
#endif