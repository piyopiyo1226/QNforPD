#pragma once
#include "mesh_struct.hpp"

//#include <vector>
namespace arcsim {

	class MeshSubset
	{
	public:
		std::vector<Node *> active_nodes;
		std::vector<Node *> support_nodes;


		std::vector<Face *> get_faces();
		std::vector<Vert *> get_verts();
		std::vector<Edge *> get_edges();
		std::vector<Node *> get_all_nodes();
		void grow(int rings);
		void rebuild_faces();
		void update_support();
		void set_flag(int flag);
		void clear_flag(int flag);

		void debug();

	private:
		//void recompute_support(std::map<Node *, int> &acc);

	};


}

