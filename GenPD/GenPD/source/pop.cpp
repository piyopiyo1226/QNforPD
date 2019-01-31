#include "pop.hpp"
namespace arcsim {



	double PopOpt::objective(const double *x) const {
		for (size_t n = 0; n < mesh.nodes.size(); n++) mesh.nodes[n]->x = x0[n] + get_subvec(x, n);
		double e = internal_energy<WS>(cloth.mesh.faces, cloth.mesh.edges);
		e += constraint_energy(cons);
		for (size_t n = 0; n < mesh.nodes.size(); n++) {
			const Node *node = mesh.nodes[n];
			e += node->m * dot(node->acceleration, node->x - x0[n]);
			e += global_index::mu * norm2(node->x - x0[n]) / 2.;
		}
		return e;
	}

}