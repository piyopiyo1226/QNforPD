#include "optim.hpp"
#include "taucs.hpp"
#include "mesh_struct.hpp"
#include "physics.hpp"

namespace arcsim {
		template <Space s>
		struct LocalOpt : public NLOpt {
			std::vector<Node *> &nodes;
			std::vector<Vec3> x0;
			std::vector<Constraint *> cons;
			mutable std::vector<Vec3> F;
			mutable SpMat<Mat3x3> J;
			const std::vector<Face *> &faces;
			const std::vector<Edge *> &edges;
			LocalOpt(std::vector<Node *> &nodes, const std::vector<Face *> &faces, const std::vector<Edge *> &edges,
				const std::vector<Constraint *> &pcons);
			void update(const double *x) const;
			void initialize(double *x) const;
			void precompute(const double *x) const;
			void finalize(const double *x) const;
			double objective(const double *x) const;
			void gradient(const double *x, double *g) const;
			bool hessian(const double *x, SpMat<double> &H) const;
		};

		template <Space s>
		LocalOpt<s>::LocalOpt( std::vector<Node *> &nodes, const std::vector<Face *> &faces, const std::vector<Edge *> &edges,
			const std::vector<Constraint *> &pcons)
			: nodes(nodes)
			, faces(faces)
			, edges(edges) {
			int nn = nodes.size();
			x0.resize(nn);
			F.resize(nn);
			J = SpMat<Mat3x3>(nn, nn);
			nvar = 3 * nn;
			for (size_t i = 0; i < nn; i++) {
				x0[i] = pos<s>(nodes[i]);
				nodes[i]->index = i;
			}
			// reduce constraints
			for (size_t n = 0; n < pcons.size(); n++) {
				for (size_t i = 0; i < nn; i++) {
					if (pcons[n]->contains(nodes[i])) {
						cons.push_back(pcons[n]);
						break;
					}
				}
			}
		}

	template <Space s>
	void LocalOpt<s>::initialize(double *x) const {
		for (int i = 0; i < nvar; i++) x[i] = 0;
	}


	template<Space s>
	void LocalOpt<s>::update(const double *x) const
	{
		for (int no = 0; no < nodes.size(); no++) pos<s>(nodes[no]) = x0[no] + get_subvec(x, no);
	}

	template <Space s>
	void LocalOpt<s>::precompute(const double *x) const
	{
		//update(x);
		for (size_t n = 0; n < nodes.size(); n++) {
			F[n] = Vec3(0);
			for (size_t jj = 0; jj < J.rows[n].entries.size(); jj++) J.rows[n].entries[jj] = Mat3x3(0);
		}
		add_internal_forces<s>(faces, edges, J, F, 0);
		add_constraint_forces(cons, J, F, 0);
	}

	template <Space s>
	double LocalOpt<s>::objective(const double *x) const {
		update(x);
		double E = internal_energy<s>(faces, edges);
		double e0 = E;
		E += constraint_energy(cons);
		double e1 = E;
		if (s == WS) {
			for (size_t n = 0; n < nodes.size(); n++) {
				const Node *node = nodes[n];
				// E += node->m * dot(node->acceleration, node->x);
			}
		}
		// cout << "internal " << e0 << " cons " << e1-e0 << " acc " << E-e1 << endl;
		return E;
	}

	template <Space s>
	void LocalOpt<s>::gradient(const double *x, double *g) const {
		for (int no = 0; no < nodes.size(); no++) {
			Vec3 f = -F[no];
			if (s == WS)
				// f += nodes[no]->m*nodes[no]->acceleration; // retain acceleration
				set_subvec(g, no, f);
			if (is_bullshit(f)) {
				std::cout << "bs alert!" << std::endl;
				//Annotation::add(nodes[no]);
				//wait_key();
				system("pause");
			}
		}
	}

	template <Space s>
	bool LocalOpt<s>::hessian(const double *x, SpMat<double> &H) const {
		double tr = 0;  // for regularization
		for (size_t i = 0; i < nodes.size(); i++) {
			const SpVec<Mat3x3> &Ji = J.rows[i];
			for (int jj = 0; jj < Ji.indices.size(); jj++) {
				int j = Ji.indices[jj];
				const Mat3x3 &Jij = Ji.entries[jj];
				set_submat(H, i, j, Jij);
			}
			tr += trace(get_submat(H, i, i));
		}
		Mat3x3 reg = Mat3x3(1e-12 * tr / nodes.size());  // regularization
		for (size_t i = 0; i < nodes.size(); i++) add_submat(H, i, i, reg);
		// TODO: should add constraint hessian, too
		return true;
	}


	template <Space s>
	void LocalOpt<s>::finalize(const double *x) const {
		update(x);
	}





	template <Space s>
	void local_opt(std::vector<Node *> &nodes, std::vector<Face *> &faces, std::vector<Edge *> &edges, const std::vector<Constraint *> &cons) {
		// mark nodes active for physics
		activate_nodes(nodes);
		line_search_newtons_method(LocalOpt<s>(nodes, faces, edges, cons), OptOptions().max_iter(10));
		deactivate_nodes(nodes);
		compute_ws_data(nodes);
		compute_ws_data(faces);
	}
	template void local_opt<PS>(std::vector<Node *> &nodes, std::vector<Face *> &faces, std::vector<Edge *> &edges,
		const std::vector<Constraint *> &cons);
	template void local_opt<WS>(std::vector<Node *> &nodes, std::vector<Face *> &faces, std::vector<Edge *> &edges,
		const std::vector<Constraint *> &cons);


	
}