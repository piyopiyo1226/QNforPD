
//------------------------------called by every---------------//
#ifndef MESH_STRUCT_HPP
#define MESH_STRUCT_HPP
#include <iostream>
#include <vector>
#include <utility>
#include "vectors.hpp"
#include "transformation.hpp"

namespace arcsim {
	//global extern
	extern int uuid_src;

	// material space (not fused at seams)
	struct Vert;
	struct Face;
	////// world space (fused)
	struct Node;
	struct Edge;
	struct Material;
	class ReferenceShape;
	struct Cloth;
	struct Mesh;
	class CollisionProxy;

	//----------------------plane_struct-----------------------//

	struct Plane
	{
		Plane() {}
		Plane(const Vec3 &x0, const Vec3 &n): 
			x0(x0), 
			n(n) {}
		Vec3 x0, n;
	};

	//----------------------vert_struct--------------//
	struct Vert
	{
		Vec3 u; //material space
		Node *node;// world space
				   //topological data
		std::vector<Face *> adjf; //adjacent face //hoshi
		int index;  //position in mesh.verts
					// derived material-space data that only changes with remeshing
					// remeshing data
		Mat3x3 sizing;// sizing matrix

		Vert()
			:node(0)//hoshi
			, index(-1) {}//hoshi	１
		explicit Vert(const Vec3 &u)
			:u(u) {}
		//void serializer(Serialize &s);

	};

	//---------------Node---------------//

	struct Node
	{
		enum NodeFlags { FlagNone = 0, FlagActive = 1, FlagMayBreak = 2, FlagResolveUni = 4, FlagResolveMax = 8 };
		int uuid;

		Mesh *mesh;

		double sep;
		int label;
		int flag;
		std::vector<arcsim::Vert *> verts;
		Vec3 y;//plastic embedding
		Vec3 x, x0, v; // position, old (collision-free) position, velocity
		bool preserve;//なんか大切らしい

		int index;
		std::vector<Edge *> adje;//adjacent edge
								  //derived world space data that change every frame
		Vec3 n;//local normal, approximate
			   //derived from material space that chnges only change with remeshing
		double a, m;//area mnass
		Mat3x3 curvature;// filtered curvature for bending fracture
						 //pop fileter data
		Vec3 acceleration;
		Node()
			: uuid(uuid_src++)
			, sep(0)
			, label(0)
			, flag(0)
			, preserve(false)
			, index(-1)
			, a(0)
			, m(0) {}
		explicit Node(const Vec3 &y, const Vec3 &x, const Vec3 &v, int label, int flag, bool preserve)
			: uuid(uuid_src++)
			, mesh(0)
			, sep(0) 
			, label(label)
			, flag(flag)
			, y(y)
			, x(x)
			, x0(x)
			, v(v)
			, preserve(preserve)
			, curvature(0) {}
		inline bool active() const { return flag & FlagActive; }

		//	void serializer(Serialize &s);

	};

	//------------------edge_struct---------//
	struct Edge
	{
		Node *n[2];//nodes
		int preserve;
		//topological data
		Face *adjf[2];
		int index;	//position in mesh.edges
					//plasticity data
		double theta_ideal, damage; //rest dihedralm damage parameer
									//constructor
		Edge()
			:index(-1)
			, theta_ideal(0)
			, damage(0)
		{
			n[0] = n[1] = NULL;
			adjf[0] = adjf[1] = NULL;
		}

		explicit Edge(Node *node0, Node *node1, double theta_ideal, int preserve)
			:preserve(preserve), theta_ideal(theta_ideal), damage(0)
		{
			n[0] = node0;
			n[1] = node1;
		}
		//	void serializer(Serialize &s);
	};

	//--------------face_struct-------------//

	struct Face {
		Vert *v[3];  // verts
		Material *material;
		int flag;
		// topological data
		Edge *adje[3];  // adjacent edges
		int index;      // position in mesh.faces
						// derived world-space data that changes every frame
		Vec3 n;  // local normal, exact
				 // derived material-space data that only changes with remeshing
		double a, m;       // area, mass
		Mat3x3 Dm, invDm;  // finite element matrix
						   // plasticity data
		Mat3x3 Sp_bend;  // plastic bending strain
		Mat3x3 Sp_str;   // plastic stretching
		Mat3x3 sigma;
		double damage;  // accumulated norm of S_plastic/S_yield
						// constructors
		Face()
			: material(0)
			, flag(0)
			, index(-1)
			, a(0)
			, m(0)
			, damage(0) {
			for (int i = 0; i < 3; i++) {
				v[i] = 0;
				adje[i] = 0;
			}
		}
		explicit Face(Vert *vert0, Vert *vert1, Vert *vert2, const Mat3x3 &ps, const Mat3x3 &pb, Material *mat,
			double damage)
			: material(mat)
			, flag(0)
			, a(0)
			, m(0)
			, Sp_bend(pb)
			, Sp_str(ps)
			, sigma(0)
			, damage(damage) {
			v[0] = vert0;
			v[1] = vert1;
			v[2] = vert2;
		}
		//void serializer(Serialize &s);
	};

	//---------------mesh_struct-----------//


	struct Mesh
	{
		ReferenceShape *ref;
		Cloth *parent;
		CollisionProxy *proxy;

		std::vector<Vert *> verts;
		std::vector<Node *> nodes;
		std::vector<Edge *> edges;
		std::vector<Face *> faces;

		void add(Vert *vert);
		void add(Node *node);
		void add(Edge *edge);
		void add(Face *face);

		void remove(Vert *vert);
		void remove(Node *node);
		void remove(Edge *edge);
		void remove(Face *face);
		Mesh()
			: ref(0)
			, parent(0)
			, proxy(0) {};

		//	void serializer(Serialize &s);

	};
	//---mesh function----------------//




	//----------------inline function--------//

	inline Vert *get_vert(const Face *face, const Node *node);
	inline Edge *get_edge(const Node *node0, const Node *node1);
	inline Vert *edge_vert(const Edge *edge, int side, int i);
	inline Vert *edge_opp_vert(const Edge *edge, int side);

	inline Vec3 derivative(double a0, double a1, double a2, double az, const Face *face);
	inline Mat3x3 derivative(const Vec3 &w0, const Vec3 &w1, const Vec3 &w2, const Vec3 &dz, const Face *face);

	//--------------------function for mesh------------//
	void delete_mesh(Mesh &mesh);
	Mesh deep_copy(Mesh &mesh);
	inline Face *adj_face(const Face *face0, int num) {
    Edge *e = face0->adje[num];
    return e->adjf[0] == face0 ? e->adjf[1] : e->adjf[0];
}
	inline Vert *get_vert(const Face *face, const Node *node) {
		if (face->v[0]->node == node) return face->v[0];
		return face->v[1]->node == node ? face->v[1] : face->v[2];
	}


	inline Edge *get_edge(const Node *n0, const Node *n1) {
		for (int e = 0; e < (int)n0->adje.size(); e++) {
			Edge *edge = n0->adje[e];
			if (edge->n[0] == n1 || edge->n[1] == n1) return edge;
		}
		return NULL;
	}

	inline Vert *edge_vert(const Edge *edge, int side, int i) {
		Face *face = (Face *)edge->adjf[side];
		if (!face) return NULL;
		for (int j = 0; j < 3; j++)
			if (face->v[j]->node == edge->n[i]) return face->v[j];
		return NULL;
	}

	inline Vert *edge_opp_vert(const Edge *edge, int side) {
		Face *face = (Face *)edge->adjf[side];
		if (!face) return NULL;
		for (int j = 0; j < 3; j++)
			if (face->v[j]->node == edge->n[side]) return face->v[j > 0 ? j - 1 : j + 2];
		return NULL;
	}

	inline Vec3 derivative(double a0, double a1, double a2, double az, const Face *face) {
		return face->invDm.t() * Vec3(a1 - a0, a2 - a0, az);
	}

	inline Mat3x3 derivative(const Vec3 &w0, const Vec3 &w1, const Vec3 &w2, const Vec3 &dz, const Face *face) {
		return Mat3x3(w1 - w0, w2 - w0, dz) * face->invDm;
	}

	template <typename Prim>
	inline const std::vector<Prim *> &get(const Mesh &mesh);
	template <typename Prim>
	inline int count_elements(const std::vector<Mesh *> &meshes);


	template <>
	inline const std::vector<Vert *> &get(const Mesh &mesh) {
		return mesh.verts;
	}
	template <>
	inline const std::vector<Node *> &get(const Mesh &mesh) {
		return mesh.nodes;
	}
	template <>
	inline const std::vector<Edge *> &get(const Mesh &mesh) {
		return mesh.edges;
	}
	template <>
	inline const std::vector<Face *> &get(const Mesh &mesh) {
		return mesh.faces;
	}

	template <typename Prim>
	inline int count_elements(const std::vector<Mesh *> &meshes) {
		int num = 0;
		for (size_t i = 0; i < meshes.size(); i++) num += get<Prim>(*meshes[i]).size();
		return num;
	}


	//-------------------------------------function
	void compute_ms_data(Face *face);
	void compute_ms_data(Node *node);
	void compute_ws_data(Node *node);
	void compute_ws_data(std::vector<Face *> &faces);
	void compute_ws_data(std::vector<Node *> &nodes);

	void set_indices(Mesh &mesh);
	void set_indices(std::vector<Mesh *> &meshes);
	void deactivate_nodes(std::vector<Node *> &nodes);


	void activate_nodes(std::vector<Node *> &nodes);

}
#endif



