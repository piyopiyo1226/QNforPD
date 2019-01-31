#include "collision.hpp"
#include <omp.h>

namespace arcsim {

	std::vector<AccelStruct * > create_accel_structs(const std::vector<Mesh * > &meshes, bool ccd)
	{
		std::vector<AccelStruct * > accs;
		for (int i = 0; i < (int)meshes.size(); i++) if (!meshes[i]->proxy)accs.push_back(new AccelStruct(*meshes[i], ccd));
		return accs;
	};

	void destroy_accel_structs(std::vector<AccelStruct *> &accs) {
		for (int a = 0; a < (int)accs.size(); a++)
			if (accs[a]) delete accs[a];
	};
	void for_overlapping_faces(BVHNode *node, float thickness, BVHCallback callback) {
		if (node->isLeaf() || !node->_active) return;
		for_overlapping_faces(node->getLeftChild(), thickness, callback);
		for_overlapping_faces(node->getRightChild(), thickness, callback);
		for_overlapping_faces(node->getLeftChild(), node->getRightChild(), thickness, callback);
	}

	void for_overlapping_faces(BVHNode *node0, BVHNode *node1, float thickness, BVHCallback callback) {
		if (!node0->_active && !node1->_active) return;
		if (!overlap(node0->_box, node1->_box, thickness)) return;
		if (node0->isLeaf() && node1->isLeaf()) {
			Face *face0 = node0->getFace(), *face1 = node1->getFace();
			callback(face0, face1);
		}
		else if (node0->isLeaf()) {
			for_overlapping_faces(node0, node1->getLeftChild(), thickness, callback);
			for_overlapping_faces(node0, node1->getRightChild(), thickness, callback);
		}
		else {
			for_overlapping_faces(node0->getLeftChild(), node1, thickness, callback);
			for_overlapping_faces(node0->getRightChild(), node1, thickness, callback);
		}
	}

	void for_overlapping_faces(const std::vector<AccelStruct *> &accs, const std::vector<AccelStruct *> &obs_accs, double thickness,
		BVHCallback callback, bool parallel, bool only_obs) {
		int nnodes = (int)ceil(sqrt(2 * omp_get_max_threads()));
		std::vector<BVHNode *> nodes = collect_upper_nodes(accs, nnodes);
		int nthreads = omp_get_max_threads();
		omp_set_num_threads(parallel ? omp_get_max_threads() : 1);
	#pragma omp parallel for
		for (int n = 0; n < (int)nodes.size(); n++) {
			if (!only_obs) {
				for_overlapping_faces(nodes[n], thickness, callback);
				for (int m = 0; m < n; m++) for_overlapping_faces(nodes[n], nodes[m], thickness, callback);
			}
			for (int o = 0; o < (int)obs_accs.size(); o++)
				if (obs_accs[o]->root) for_overlapping_faces(nodes[n], obs_accs[o]->root, thickness, callback);
		}
		omp_set_num_threads(nthreads);
	}


	std::vector<BVHNode *> collect_upper_nodes(const std::vector<AccelStruct *> &accs, int nnodes) {
		std::vector<BVHNode *> nodes;
		for (int a = 0; a < (int)accs.size(); a++)
			if (accs[a]->root) nodes.push_back(accs[a]->root);
		while ((int)nodes.size() < nnodes) {
			std::vector<BVHNode *> children;
			for (int n = 0; n < (int)nodes.size(); n++)
				if (nodes[n]->isLeaf())
					children.push_back(nodes[n]);
				else {
					children.push_back(nodes[n]->_left);
					children.push_back(nodes[n]->_right);
				}
				if (children.size() == nodes.size()) break;
				nodes = children;
		}
		return nodes;
	}



	void collect_leaves(BVHNode *node, std::map<const Face *, BVHNode *> &leaves) {
		if (node->isLeaf()) {
			leaves[node->getFace()] = node;
		}
		else {
			collect_leaves(node->getLeftChild(), leaves);
			collect_leaves(node->getRightChild(), leaves);
		}
	}


	//momo
	namespace global_index{
		const std::vector<Mesh *> *meshes, *obs_meshes;
	}
}