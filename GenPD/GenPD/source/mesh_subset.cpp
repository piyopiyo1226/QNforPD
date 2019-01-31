#include "mesh_subset.hpp"
#include "util.hpp"
namespace arcsim {

	std::vector<Face *> MeshSubset::get_faces() {
		std::vector<Face *> faces;
		for (size_t i = 0; i < active_nodes.size(); i++) {
			std::vector<Vert *> &verts = active_nodes[i]->verts;
			for (size_t j = 0; j < verts.size(); j++) {
				std::vector<Face *> &f = verts[j]->adjf;
				for (size_t k = 0; k < f.size(); k++) include(f[k], faces);
			}
		}
		return faces;
	}
}