#include "mesh_struct.hpp"
#include "Sparce.hpp"
#include "geometry.hpp"

#include "arcsim_constraint.hpp"//constraint
#include <vector>
#include <iostream>
namespace arcsim {

	template <Space s>
	void local_opt(std::vector<Node *> &nodes, std::vector<Face *> &faces, std::vector<Edge *> &edges,const std::vector<Constraint *> &cons);

	struct NLOpt {  // nonlinear optimization problem
					// minimize objective
		int nvar;//‚È‚É‚©‚È
		virtual void initialize(double *x) const = 0;
		virtual double objective(const double *x) const = 0;
		virtual void precompute(const double *x) const {}
		virtual void gradient(const double *x, double *g) const = 0;
		virtual bool hessian(const double *x, SpMat<double> &H) const {
			return false;  // should return true if implemented
		};
		virtual void finalize(const double *x) const = 0;
	};





	
	struct OptOptions
	{
		int _max_iter;
		double _eps_x, _eps_f, _eps_g;
		OptOptions()
			: _max_iter(100)
			, _eps_x(1e-6)
			, _eps_f(1e-12)
			, _eps_g(1e-6) {}
		// Named parameter idiom
		// http://www.parashift.com/c++-faq-lite/named-parameter-idiom.html
		OptOptions &max_iter(int n) {
			_max_iter = n;
			return *this;
		}
		OptOptions &eps_x(double e) {
			_eps_x = e;
			return *this;
		}
		OptOptions &eps_f(double e) {
			_eps_f = e;
			return *this;
		}
		OptOptions &eps_g(double e) {
			_eps_g = e;
			return *this;
		}
		int max_iter() { return _max_iter; }
		double eps_x() { return _eps_x; }
		double eps_f() { return _eps_f; }
		double eps_g() { return _eps_g; }
	};


	void line_search_newtons_method(const NLOpt &problem, OptOptions opts = OptOptions(), bool verbose = false);


	inline void set_subvec(double *x, int i, const Vec3 &xi) {
		for (int j = 0; j < 3; j++) x[i * 3 + j] = xi[j];
	}

	template <int n>
	void set_subvec(double *x, int i, const Vec<n> &xi) {
		for (int j = 0; j < n; j++) x[i * n + j] = xi[j];
	}

	inline void set_submat(SpMat<double> &A, int i, int j, const Mat3x3 &Aij) {
		for (int ii = 0; ii < 3; ii++)
			for (int jj = 0; jj < 3; jj++) A(i * 3 + ii, j * 3 + jj) = Aij(ii, jj);
	}
	inline Mat3x3 get_submat(SpMat<double> &A, int i, int j) {
		Mat3x3 Aij;
		for (int ii = 0; ii < 3; ii++)
			for (int jj = 0; jj < 3; jj++) Aij(ii, jj) = A(i * 3 + ii, j * 3 + jj);
		return Aij;
	}
	inline void add_submat(SpMat<double> &A, int i, int j, const Mat3x3 &Aij) {
		for (int ii = 0; ii < 3; ii++)
			for (int jj = 0; jj < 3; jj++) A(i * 3 + ii, j * 3 + jj) += Aij(ii, jj);
	}
}

//

