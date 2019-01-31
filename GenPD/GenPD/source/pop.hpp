#ifndef POPFILTER_HPP
#define POPFILTER_HPP

#include "sub_struct.hpp"
#include "arcsim_constraint.hpp"
#include "optim.hpp"
#include "sparce.hpp"
#include "physics.hpp"

namespace arcsim {
	namespace global_index {
		static double mu;
	}
	struct PopOpt : public NLOpt {
		// we want F(x) = m a0, but F(x) = -grad E(x)
		// so grad E(x) + m a0 = 0
		// let's minimize E(x) + m a0 . (x - x0)
		// add spring to x0: minimize E(x) + m a0 . (x - x0) + mu (x - x0)^2/2
		// gradient: -F(x) + m a0 + mu (x - x0)
		// hessian: -J(x) + mu
		Cloth &cloth;
		Mesh &mesh;
		const std::vector<Constraint *> &cons;
		std::vector<Vec3> x0, a0;
		mutable std::vector<Vec3> f;
		mutable SpMat<Mat3x3> J;
		PopOpt(Cloth &cloth, const std::vector<Constraint *> &cons)
			: cloth(cloth)
			, mesh(cloth.mesh)
			, cons(cons) {
			int nn = mesh.nodes.size();
			nvar = nn * 3;
			x0.resize(nn);
			a0.resize(nn);
			for (int n = 0; n < nn; n++) {
				const Node *node = mesh.nodes[n];
				x0[n] = node->x;
				a0[n] = node->acceleration;
			}
			f.resize(nn);
			J = SpMat<Mat3x3>(nn, nn);
		}
		virtual void initialize(double *x) const;
		virtual void precompute(const double *x) const;
		virtual double objective(const double *x) const;
		virtual void gradient(const double *x, double *g) const;
		virtual bool hessian(const double *x, SpMat<double> &H) const;
		virtual void finalize(const double *x) const;
	};
	void apply_pop_filter(Cloth &cloth, const std::vector<Constraint *> &cons, double regularization = 1e3);
}
#endif
