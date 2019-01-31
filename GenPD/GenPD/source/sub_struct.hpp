#ifndef SUB_HPP
#define SUB_HPP
#include "mesh_struct.hpp"
//#include "arcsim_constraint.hpp"
//#include "magic.hpp"
#include "collision.hpp"
#include "dde.hpp"

namespace arcsim {

	struct Material
	{
		double density;//area density
		StretchingSamples dde_stretching;
		BendingData dde_bending;
		double damping;//stiffness-proportional damping coefficient
		double strain_min, strain_max;//strain max min
		double yield_curv, weakening; //plasticity parameters
		double yield_stretch, plastic_flow, plastic_limit;
		bool use_dde;//use dde material files
		double thickness;
		double alt_stretching, alt_bending, alt_poisson;  // alternative material model
		double toughness, fracture_bend_thickness;
	};

	struct Remeshing
	{
		double refine_angle, refine_compression, refine_velocity;
		double size_min, size_max, size_uniform;
		double aspect_min;
		double refine_fracture;
	};

	struct Cloth {//
		Mesh mesh;
		std::vector<Material  *> materials;
		Remeshing remeshing;
	};

	void compute_material(Material &mat, double Y);
	//---------------------------------------------------------------------//

	

	struct Constraint;

	
	// Proxy for obstacle collisions. Helpful for fast-moving
	class CollisionProxy {
	public:
		virtual ~CollisionProxy() {};
		virtual CollisionProxy *clone(Mesh &mesh) = 0;
		virtual void update(Mesh &mesh) = 0;
		virtual Constraint *constraint(const Node *node) = 0;
	};


	//--------------------proxy hpp---------------------------//

	class FloorProxy : public CollisionProxy {
	public:
		FloorProxy(Mesh &mesh);

		Constraint *constraint(const Node *node);
		CollisionProxy *clone(Mesh &mesh);
		void update(Mesh &mesh);

	private:
		Node center;
	};

}
//---------------------------sim_struct----------
#endif