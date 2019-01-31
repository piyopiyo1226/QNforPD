#pragma once

#include "sub_struct.hpp"
#include "mesh_subset.hpp"

#include <map>
#include <algorithm>
#include <vector>

namespace arcsim {


	// Helpers---for remesh.hpp--------------------------------
	template <class T>
	static void delete_all(const std::vector<T> &a) {
		for (size_t i = 0; i < a.size(); i++) delete a[i];
	}
	template <class T>
	static void remove_all(const std::vector<T> &a, Mesh &m) {
		for (size_t i = 0; i < a.size(); i++) m.remove(a[i]);
	}
	template <class T>
	static void add_all(const std::vector<T> &a, Mesh &m) {
		for (size_t i = 0; i < a.size(); i++) m.add(a[i]);
	}
	template <class T>
	static void include_all(const std::vector<T> &a, std::vector<T> &b) {
		for (size_t i = 0; i < a.size(); i++) include(a[i], b);
	}
	template <class T>
	static void exclude_all(const std::vector<T> &a, std::vector<T> &b) {
		for (size_t i = 0; i < a.size(); i++) exclude(a[i], b);
	}

	struct RemeshOp {
		std::vector<Vert *> added_verts, removed_verts;
		std::vector<Node *> added_nodes, removed_nodes;
		std::vector<Edge *> added_edges, removed_edges;
		std::vector<Face *> added_faces, removed_faces;

		bool empty() { return added_faces.empty() && removed_faces.empty(); }
		RemeshOp inverse() const;
		void apply(Mesh &mesh) const;
		void done() const;  // frees removed data
		void cancel();      // aborts operations, resets RemeshOp

		void set_null(std::vector<Edge *> &v);
		void update(std::vector<Node *> &v);
		void update(std::vector<Face *> &v);
		void update(std::vector<Edge *> &v);
	};





	void dynamic_remesh(Mesh &mesh, const std::map<Node *, Plane> &planes);
	void delete_spaced_out(Mesh &mesh);

	Mat3x3 compute_face_sizing(Remeshing &remeshing, const Face *face, const std::map<Node *, Plane> &planes,
		bool debug = false);
	Mat2x2 compression_metric(const Face *face, const Mat3x3 &S2, const Mat3x2 &UV, double c);

	//---------------------------for function------------------------------//

	void flip_edges(MeshSubset *subset, std::vector<Face *> &active_faces, std::vector<Edge *> *update_edges, std::vector<Face *> *update_faces);
	//	bool try_move_node(Node *node, Edge *edge, double d);
	RemeshOp split_edge(Edge *edge, double d);
	RemeshOp collapse_edge(Edge *edge, int which);  // which end to delete
	RemeshOp flip_edge(Edge *edge);

	void  create_vert_sizing(std::vector<Vert *> &verts, const std::map<Node *, Plane> &planes);
	bool split_worst_edge(MeshSubset *subset, const std::vector<Edge *> &edges);
	bool improve_some_face(MeshSubset *subset, std::vector<Face *> &active);

	Mat3x3 obstacle_metric(const Face *face, const std::map<Node *, Plane> &planes);
	Mat2x2 fracture_metric(Remeshing &remeshing, const Face *face);


	void embedding_from_plasticity(const std::vector<Face *> &fs);
	void plasticity_from_embedding(const std::vector<Face *> faces);

	void local_pop_filter(const std::vector<Face * > &fs);

	std::vector<Edge *> find_edges_to_flip(std::vector<Face *> &active_faces);
	bool should_flip(const Edge *edge);

	std::vector<Edge *> find_bad_edges(const std::vector<Edge *> &edges);
	Vert *adjacent_vert(const Node *node, const Vert *vert);

	bool has_labeled_edges(const Node *node);
	RemeshOp try_edge_collapse(Edge *edge, int which);
	bool can_collapse(Remeshing &remeshing, const Edge *edge, int i);
		


		void optimize_node(Node *node);
	//	Mat3x3 stretch_plasticity_from_embedding(const Face *face);
		void recompute_Sp_bend(Face *face);
		

}

