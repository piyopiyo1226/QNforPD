#include "sub_struct.hpp"
#include "arcsim_simulation.hpp"
#include "magic.hpp"

namespace arcsim{

void compute_material(Material &mat, double Y) {
	double A = Y / (1.0 - sq(mat.alt_poisson));
	mat.alt_stretching = A * mat.thickness;
	mat.alt_bending = A / 12.0 * mat.thickness * sq(mat.thickness);
	mat.toughness *= mat.thickness;
}


//-------------------for proxy-----//


FloorProxy::FloorProxy(Mesh &mesh) { update(mesh); }

CollisionProxy *FloorProxy::clone(Mesh &mesh) { return new FloorProxy(mesh); }
//
Constraint *FloorProxy::constraint(const Node *node) {
	if (!is_free(node) || node->x[1] - center.x[1] > global_index::magic.repulsion_thickness) return 0;

	IneqCon *con = new IneqCon;
	con->nodes[0] = (Node *)node;
	con->nodes[1] = &center;
	con->nodes[2] = 0;
	con->nodes[3] = 0;
	for (int n = 0; n < 4; n++) con->free[n] = n == 0;
	con->w[0] = 1;
	con->w[1] = -1;
	con->w[2] = 0;
	con->w[3] = 0;

	con->stiff = global_index::magic.collision_stiffness * node->a;
	con->n = Vec3(0, 1, 0);

	con->mu = sim.obs_friction;
	return con;
}

void FloorProxy::update(Mesh &mesh) {
	double y = 0;// infinity;//hosi
	for (size_t n = 0; n < mesh.nodes.size(); n++)
		if (mesh.nodes[n]->x[1] < y) y = mesh.nodes[n]->x[1];
	center.x = Vec3(0, y, 0);
}
}
