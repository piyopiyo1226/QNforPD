#include "optim.hpp"
#include "taucs.hpp"
namespace arcsim
{
	//‚È‚ºstatic
	namespace global_index {
		static bool verbose;
	}


//	static double dot(const std::vector<double> &x, const std::vector<double> &y);
	static double dot(const std::vector<double> &x, const std::vector<double> &y) {
		int n = x.size();
		double d = 0;
		for (int i = 0; i < n; i++) d += x[i] * y[i];
		return d;
	}
	static double norm(const std::vector<double> &x) { return sqrt(dot(x, x)); }
	static void scalar_mult(std::vector<double> &v, double a, const std::vector<double> &x) {
		int n = x.size();
		v.resize(n);
		for (int i = 0; i < n; i++) v[i] = a * x[i];
	}


	static void add(std::vector<double> &v, double a, const std::vector<double> &x, double b, const std::vector<double> &y) {
		int n = x.size();
		v.resize(n);
		for (int i = 0; i < n; i++) v[i] = a * x[i] + b * y[i];
	}

	inline double cb(double x) { return x * x * x; }


	double line_search(const std::vector<double> &x0, const std::vector<double> &p, const NLOpt &problem, double f0,
		const std::vector<double> &g) {
		double c = 1e-3;  // sufficient decrease parameter
		double a = 1;
		int n = problem.nvar;
		std::vector<double> x(n);
		double g0 = dot(g, p);
		if (global_index::verbose) REPORT(g0);
		if (abs(g0) < 1e-12) return 0;
		double a_prev = 0;
		double f_prev = f0;
		while (true) {
			add(x, 1, x0, a, p);
			// problem.precompute(&x[0]);
			double f = problem.objective(&x[0]);
			if (global_index::verbose) {
				REPORT(a);
				REPORT(f);
			}
			if (f <= f0 + c * a * g0) break;
			double a_next;
			if (a_prev == 0)
				a_next = -g0 * sq(a) / (2 * (f - f0 - g0 * a));  // minimize quadratic fit
			else {
				// minimize cubic fit to f0, g0, f_prev, f
				Vec2 b = Mat2x2(Vec2(sq(a_prev), -cb(a_prev)), Vec2(-sq(a), cb(a))) *
					Vec2(f - f0 - g0 * a, f_prev - f0 - g0 * a_prev) / (sq(a) * sq(a_prev) * (a - a_prev));
				double a_sol[2];
				solve_quadratic(3 * b[0], 2 * b[1], g0, a_sol);
				a_next = (a_sol[0] > 0) ? a_sol[0] : a_sol[1];
			}
			if (a_next < a * 0.1 || a_next > a * 0.9) a_next = a / 2;
			a_prev = a;
			f_prev = f;
			a = a_next;
		}
		return a;
	}

	void line_search_newtons_method(const NLOpt &problem, OptOptions opt, bool verbose) {//check for optim NLOpt and OptOptions
		global_index::verbose = verbose;//hoshi
		int n = problem.nvar;
		std::vector<double> x(n), g(n);
		SpMat<double> H(n, n);
		problem.initialize(&x[0]);
		double f_old = infinity;
		int iter;
		for (iter = 0; iter < opt.max_iter(); iter++) {
			problem.precompute(&x[0]);
			double f = problem.objective(&x[0]);
			if (verbose) REPORT(f);
			if (f_old - f < opt.eps_f()) break;
			f_old = f;
			problem.gradient(&x[0], &g[0]);
			if (verbose) REPORT(norm(g));//vectors.hpp
			if (norm(g) < opt.eps_g()) break;
			if (!problem.hessian(&x[0], H)) {
				std::cerr << "Can't run Newton's method if Hessian of objective "
					<< "is not available!" << std::endl;
				exit(1);
			}
			std::vector<double> p = taucs_linear_solve(H, g);//rewrite the taucs.hpp folder
			if (verbose) REPORT(norm(p));
			scalar_mult(p, -1, p);
			double a = line_search(x, p, problem, f, g);//4-12-1
			add(x, 1, x, a, p);
			if (a * norm(p) < opt.eps_x()) break;
		}
		if (verbose) REPORT(iter);
		problem.finalize(&x[0]);
	}

}
