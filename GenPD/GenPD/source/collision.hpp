#pragma once
//-----------------in collision hpp----------------//
#include "bvh.hpp"

#include <map>

namespace arcsim {
	typedef DeformBVHNode BVHNode;
	typedef DeformBVHTree BVHTree;

	// callback must be safe to parallelize via OpenMP
	typedef void(*BVHCallback)(const Face *face0, const Face *face1);
	void collect_leaves(BVHNode *node, std::map<const Face *, BVHNode *> &leaves);


	struct AccelStruct
	{
		BVHTree tree;
		BVHNode *root;
		std::map<const Face*, BVHNode *> leaves;
		//	AccelStruct(const Mesh &mesh, bool ccd);

		AccelStruct(const Mesh &mesh, bool ccd) :tree((Mesh &)mesh, ccd), root(tree._root), leaves()
		{
			if (root) { collect_leaves(root, leaves); };
		}
	};

	

	void destroy_accel_structs(std::vector<AccelStruct *> &accs);
	


	//-------------
	std::vector<AccelStruct * > create_accel_structs(const std::vector<Mesh * > &meshes, bool ccd);
	void for_overlapping_faces(const std::vector<AccelStruct *> &accs, const std::vector<AccelStruct *> &obs_accs, double thickness, BVHCallback callback, bool parallel, bool only_obs);
	void for_overlapping_faces(BVHNode *node, float thickness, BVHCallback callback);
	void for_overlapping_faces(BVHNode *node0, BVHNode *node1, float thickness, BVHCallback callback);
	std::vector<BVHNode *> collect_upper_nodes(const std::vector<AccelStruct *> &accs, int nnodes);


//collisiont util
	template <typename T>
	inline bool is_free(const T *p);
	
	template<>
	inline bool is_free <Node>(const Node *p)
	{
		return p->verts[0]->adjf[0]->material;
	}
	template <>
	inline bool is_free<Vert>(const Vert *p) {
		return p->adjf[0]->material;
	}
	template <>
	inline bool is_free<Edge>(const Edge *p) {
		return p->n[0]->verts[0]->adjf[0]->material;
	}
	template <>
	inline bool is_free<Face>(const Face *p) {
		return p->material;
	}

	namespace global_index {
		extern const std::vector<Mesh *> *meshes;
	}
}


