#include "proximity.hpp"

namespace arcsim
{
	namespace global_index {
		static std::vector<Min<Face *> > node_prox[2];
		static std::vector<Min<Edge *> > edge_prox[2];
		static std::vector<Min<Node *> > face_prox[2];
		static std::vector<Min<Node *> > edge_node_prox;
		static std::vector<Min<Edge *> > node_edge_prox;
	}

	double area_cached(const Node *node);
	double area_cached(const Edge *edge);
	double area_cached(const Face *face);

	static inline bool has_node(const Face *f, const Node *n) {
		return f && (f->v[0]->node == n || f->v[1]->node == n || f->v[2]->node == n);
	}

	static inline Vec3 get_outwards_normal(const Edge *e) {
		Vec3 p0 = e->n[0]->x, p1 = e->n[1]->x, ne = normalize(p1 - p0);
		Vec3 opp = p0 - edge_opp_vert(e, e->adjf[0] ? 0 : 1)->node->x;
		Vec3 n = cross(e->adjf[0] ? e->adjf[0]->n : e->adjf[1]->n, ne);
		if (dot(n, opp) < 0) n = -n;
		return n;
	}


	std::vector<Constraint *> proximity_constraints(std::vector<Mesh *> &meshes, const std::vector<Mesh *> &obs_meshes, double mu,
		double mu_obs, bool proxy_only, bool obs_only) {
		global_index::meshes = &meshes;
		const double dmin = 2 * global_index::magic.repulsion_thickness;
		std::vector<Constraint *> cons;

		if (proxy_only || obs_only) {
			for (size_t m = 0; m < meshes.size(); m++) {
				for (size_t i = 0; i < obs_meshes.size(); i++) {
					if (obs_meshes[i]->proxy)
						make_proxy_constraints(*meshes[m], *(CollisionProxy *)(obs_meshes[i]->proxy), cons);//8-1-1
				}
			}
			if (proxy_only) return cons;
		}

		std::vector<AccelStruct *> accs = create_accel_structs(meshes, false),
			obs_accs = create_accel_structs(obs_meshes, false);

		set_indices(meshes);//8-1-2
		int nn = 0, ne = 0, nf = 0;
		for (size_t m = 0; m < meshes.size(); m++) {
			nn += meshes[m]->nodes.size();
			ne += meshes[m]->edges.size();
			nf += meshes[m]->faces.size();
		}
		for (int i = 0; i < 2; i++) {
			global_index::node_prox[i].assign(nn, Min<Face *>());
			global_index::edge_prox[i].assign(ne, Min<Edge *>());
			global_index::face_prox[i].assign(nf, Min<Node *>());
		}
		global_index::edge_node_prox.assign(ne, Min<Node *>());
		global_index::node_edge_prox.assign(nn, Min<Edge *>());

		for_overlapping_faces(accs, obs_accs, dmin, find_proximities, true, obs_only);//8-1-3,8-1-4

		for (size_t m = 0; m < meshes.size(); m++) {
			Mesh &mesh = *meshes[m];

			for (size_t n = 0; n < mesh.nodes.size(); n++) {
				int idx = mesh.nodes[n]->index;
				for (int i = 0; i < 2; i++) {
					Min<Face *> &m = global_index::node_prox[i][idx];
					if (m.key < dmin) cons.push_back(make_constraint(mesh.nodes[n], m.val, mu, mu_obs));//8-2-1
				}
				Min<Edge *> &me = global_index::node_edge_prox[idx];
				if (me.key < dmin) cons.push_back(make_constraint(me.val, mesh.nodes[n], mu, mu_obs));
			}
			for (size_t e = 0; e < mesh.edges.size(); e++) {
				int idx = mesh.edges[e]->index;
				for (int i = 0; i < 2; i++) {
					Min<Edge *> &m = global_index::edge_prox[i][idx];
					if (m.key < dmin) cons.push_back(make_constraint(mesh.edges[e], m.val, mu, mu_obs));
				}
				Min<Node *> &me = global_index::edge_node_prox[idx];
				if (me.key < dmin) cons.push_back(make_constraint(mesh.edges[e], me.val, mu, mu_obs));
			}
			for (size_t f = 0; f < mesh.faces.size(); f++) {
				int idx = mesh.faces[f]->index;
				for (int i = 0; i < 2; i++) {
					Min<Node *> &m = global_index::face_prox[i][idx];
					if (m.key < dmin) cons.push_back(make_constraint(m.val, mesh.faces[f], mu, mu_obs));
				}
			}

			for (size_t i = 0; i < obs_meshes.size(); i++) {
				if (obs_meshes[i]->proxy) make_proxy_constraints(mesh, *(CollisionProxy *)(obs_meshes[i]->proxy), cons);
			}
		}

		//if (consistency_check) {//
		//	std::cout << "> proximity2: " << std::endl;
		//	//test_state(cons, "/tmp/cs");
		//}

		destroy_accel_structs(accs);
		destroy_accel_structs(obs_accs);
		return cons;
	}

	void make_proxy_constraints(Mesh &mesh, CollisionProxy &proxy, std::vector<Constraint *> &cons) {//8-1-1
		for (size_t i = 0; i < mesh.nodes.size(); i++) {
			Constraint *c = proxy.constraint(mesh.nodes[i]);
			if (c) cons.push_back(c);
		}
	}

	void find_proximities(const Face *face0, const Face *face1) {//8-1-4
		for (int v = 0; v < 3; v++) add_proximity(face0->v[v]->node, face1);//8-1-8
		for (int v = 0; v < 3; v++) add_proximity(face1->v[v]->node, face0);
		for (int e0 = 0; e0 < 3; e0++)
			for (int e1 = 0; e1 < 3; e1++) {
				add_proximity(face0->adje[e0], face1->adje[e1]);//8-1-7
				add_proximity(face0->v[e0]->node, face1->adje[e1]);
				add_proximity(face1->v[e0]->node, face0->adje[e1]);
			}
	}
	bool in_wedge(double w, const Edge *edge0, const Edge *edge1) {//8-1-6
		Vec3 x = (1 - w) * edge0->n[0]->x + w * edge0->n[1]->x;
		bool in = true;
		for (int s = 0; s < 2; s++) {
			const Face *face = edge1->adjf[s];
			if (!face) continue;
			const Node *node0 = edge1->n[s], *node1 = edge1->n[1 - s];
			Vec3 e = node1->x - node0->x, n = face->n, r = x - node0->x;
			in &= right_handed(e, n, r);
		}
		return in;
	}

	void add_proximity(const Node *node, const Face *face) {
		if (has_node(face, node)) return;
		Vec3 n;
		double w[4];
		double d = signed_vf_distance(node->x, face->v[0]->node->x, face->v[1]->node->x, face->v[2]->node->x, &n, w);
		d = abs(d);
		bool inside = (min(-w[1], -w[2], -w[3]) >= -1e-6);
		if (!inside) return;
		if (is_free(node)) {
			int side = dot(n, node->n) >= 0 ? 0 : 1;
			global_index::node_prox[side][node->index].add(d, (Face *)face);
		}
		if (is_free(face)) {
			int side = dot(-n, face->n) >= 0 ? 0 : 1;
			global_index::face_prox[side][face->index].add(d, (Node *)node);
		}
	}
	void add_proximity(const Edge *edge0, const Edge *edge1) {
		if (edge0->n[0] == edge1->n[0] || edge0->n[0] == edge1->n[1] || edge0->n[1] == edge1->n[0] ||
			edge0->n[1] == edge1->n[1])
			return;
		Vec3 n;
		double w[4];
		double d = signed_ee_distance(edge0->n[0]->x, edge0->n[1]->x, edge1->n[0]->x, edge1->n[1]->x, &n, w);
		d = abs(d);
		bool inside =
			(min(w[0], w[1], -w[2], -w[3]) >= -1e-6 && in_wedge(w[1], edge0, edge1) && in_wedge(-w[3], edge1, edge0));
		if (!inside) return;
		if (is_free(edge0)) {
			Vec3 edge0n = edge0->n[0]->n + edge0->n[1]->n;
			int side = dot(n, edge0n) >= 0 ? 0 : 1;
			global_index::edge_prox[side][edge0->index].add(d, (Edge *)edge1);
		}
		if (is_free(edge1)) {
			Vec3 edge1n = edge1->n[0]->n + edge1->n[1]->n;
			int side = dot(-n, edge1n) >= 0 ? 0 : 1;
			global_index::edge_prox[side][edge1->index].add(d, (Edge *)edge0);
		}
	}

	void add_proximity(Node *node, Edge *edge) {
		if (node == edge->n[0] || node == edge->n[1] || !is_seam_or_boundary(node) || !is_seam_or_boundary(edge) ||
			(edge->adjf[0] && has_node(edge->adjf[0], node)) || (edge->adjf[1] && has_node(edge->adjf[1], node)))
			return;

		Vec3 normal = get_outwards_normal(edge);
		Vec3 p0 = edge->n[0]->x, p1 = edge->n[1]->x, x = node->x;
		if (dot(x - p0, normal) < 0) return;

		double d = dot(x - p0, p1 - p0) / norm2(p1 - p0);
		if (d < 0 || d > 1) return;

		double w0 = 1.0 - d, w1 = d;
		double dist = norm(w0 * p0 + w1 * p1 - x);
		if (is_free(node)) global_index::node_edge_prox[node->index].add(dist, edge);
		if (is_free(edge)) global_index::edge_node_prox[edge->index].add(dist, node);
	}


	Constraint *make_constraint(const Edge *edge0, const Edge *edge1, double mu, double mu_obs) {
		IneqCon *con = new IneqCon;
		con->nodes[0] = (Node *)edge0->n[0];
		con->nodes[1] = (Node *)edge0->n[1];
		con->nodes[2] = (Node *)edge1->n[0];
		con->nodes[3] = (Node *)edge1->n[1];
		for (int n = 0; n < 4; n++) con->free[n] = is_free(con->nodes[n]);
		double a = std::min(area_cached(edge0), area_cached(edge1));
		con->stiff = global_index::magic.collision_stiffness * a;
		double d =
			signed_ee_distance(con->nodes[0]->x, con->nodes[1]->x, con->nodes[2]->x, con->nodes[3]->x, &con->n, con->w);
		if (d < 0) con->n = -con->n;
		con->mu = (!is_free(edge0) || !is_free(edge1)) ? mu_obs : mu;
		return con;
	}

	Constraint *make_constraint(const Edge *edge, const Node *node, double mu, double mu_obs) {
		IneqCon *con = new IneqCon;
		con->nodes[0] = (Node *)node;
		con->nodes[1] = (Node *)edge->n[0];
		con->nodes[2] = (Node *)edge->n[1];
		con->nodes[3] = 0;
		for (int n = 0; n < 4; n++) con->free[n] = con->nodes[n] ? is_free(con->nodes[n]) : false;

		double a = (std::min)(area_cached(edge), area_cached(node));
		con->stiff = global_index::magic.collision_stiffness * a;
		con->n = get_outwards_normal(edge);
		double d = signed_ve_distance(con->nodes[0]->x, con->nodes[1]->x, con->nodes[2]->x, &con->n, con->w);

		if (fabs(d) > 2.0 * global_index::magic.repulsion_thickness) return 0;

		con->mu = (!is_free(node) || !is_free(edge)) ? mu_obs : mu;
		return con;
	}

	Constraint *make_constraint(const Node *node, const Face *face, double mu, double mu_obs) {//8-2-1
		IneqCon *con = new IneqCon;
		con->nodes[0] = (Node *)node;
		con->nodes[1] = (Node *)face->v[0]->node;
		con->nodes[2] = (Node *)face->v[1]->node;
		con->nodes[3] = (Node *)face->v[2]->node;
		for (int n = 0; n < 4; n++) con->free[n] = is_free(con->nodes[n]);
		double a = std::min(area_cached(node), area_cached(face));
		con->stiff = global_index::magic.collision_stiffness * a;
		double d =
			signed_vf_distance(con->nodes[0]->x, con->nodes[1]->x, con->nodes[2]->x, con->nodes[3]->x, &con->n, con->w);
		if (d < 0) con->n = -con->n;
		con->mu = (!is_free(node) || !is_free(face)) ? mu_obs : mu;
		return con;
	}
	

	double area_cached(const Node *node) {
		if (is_free(node)) return node->a;
		double a = 0;
		for (int v = 0; v < (int)node->verts.size(); v++)
			for (int f = 0; f < (int)node->verts[v]->adjf.size(); f++) a += area(node->verts[v]->adjf[f]) / 3;
		return a;
	}
	double area_cached(const Edge *edge) {
		double a = 0;
		if (edge->adjf[0]) a += area(edge->adjf[0]) / 3;
		if (edge->adjf[1]) a += area(edge->adjf[1]) / 3;
		return a;
	}

	double area_cached(const Face *face) {
		if (is_free(face)) return face->a;
		const Vec3 &x0 = face->v[0]->node->x, &x1 = face->v[1]->node->x, &x2 = face->v[2]->node->x;
		return norm(cross(x1 - x0, x2 - x0)) / 2;
	}


}