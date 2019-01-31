#pragma once
#include "taucs.hpp"
#include "util.hpp"
#include "alglib/solvers.h"
#include <cstdlib>
//using namespace alglib;



//-------------------tausc_ssc_matrix
	extern "C" {
		#include "taucs.h"
		int taucs_linsolve(taucs_ccs_matrix *A,   // input matrix
			void **factorization,  // an approximate inverse
			int nrhs,              // number of right-hand sides
			void *X,               // unknowns
			void *B,               // right-hand sides
			char *options[],       // options (what to do and how)
			void *arguments[]);    // option arguments
	}

	std::ostream &operator<<(std::ostream &out, taucs_ccs_matrix *A) {
		out << "n: " << A->n << std::endl;
		out << "m: " << A->m << std::endl;
		out << "flags: " << A->flags << std::endl;
		out << "colptr: ";
		for (int i = 0; i <= A->n; i++) out << (i == 0 ? "" : ", ") << A->colptr[i];
		out << std::endl;
		out << "rowind: ";
		for (int j = 0; j <= A->colptr[A->n]; j++) out << (j == 0 ? "" : ", ") << A->rowind[j];
		out << std::endl;
		out << "values.d: ";
		for (int j = 0; j <= A->colptr[A->n]; j++) out << (j == 0 ? "" : ", ") << A->values.d[j];
		out << std::endl;
		return out;
	}

	//-------------------------------------------1.2------------------------------------------------
	taucs_ccs_matrix *sparse_to_taucs(const arcsim::SpMat<double> &As) {
		// assumption: A is square and symmetric
		int n = As.n;
		int nnz = 0;
		for (int i = 0; i < n; i++) {
			for (int k = 0; k < (int)As.rows[i].indices.size(); k++) {
				int j = As.rows[i].indices[k];
				if (j < i) continue;
				nnz++;
			}
		}
		taucs_ccs_matrix *At = taucs_ccs_create(n, n, nnz, TAUCS_DOUBLE | TAUCS_SYMMETRIC | TAUCS_LOWER);
		int pos = 0;
		for (int i = 0; i < n; i++) {
			At->colptr[i] = pos;
			for (int k = 0; k < (int)As.rows[i].indices.size(); k++) {
				int j = As.rows[i].indices[k];
				if (j < i) continue;
				At->rowind[pos] = j;
				At->values.d[pos] = As.rows[i].entries[k];
				pos++;
			}
		}
		At->colptr[n] = pos;
		return At;
	}
	template <int m>
	taucs_ccs_matrix *sparse_to_taucs(const arcsim::SpMat<arcsim::Mat<m, m> > &As) {
		// assumption: A is square and symmetric
		int n = As.n;
		int nnz = 0;
		for (int i = 0; i < n; i++) {
			for (int jj = 0; jj < (int)As.rows[i].indices.size(); jj++) {
				int j = As.rows[i].indices[jj];
				if (j < i) continue;
				nnz += (j == i) ? m * (m + 1) / 2 : m * m;
			}
		}
		taucs_ccs_matrix *At = taucs_ccs_create(n * m, n * m, nnz, TAUCS_DOUBLE | TAUCS_SYMMETRIC | TAUCS_LOWER);
		int pos = 0;
		for (int i = 0; i < n; i++) {
			for (int k = 0; k < m; k++) {
				At->colptr[i * m + k] = pos;
				for (int jj = 0; jj < (int)As.rows[i].indices.size(); jj++) {
					int j = As.rows[i].indices[jj];
					if (j < i) continue;
					const arcsim::Mat<m, m> &Aij = As.rows[i].entries[jj];
					for (int l = (i == j) ? k : 0; l < m; l++) {
						At->rowind[pos] = j * m + l;
						At->values.d[pos] = Aij(k, l);
						pos++;
					}
				}
			}
		}
		At->colptr[n * m] = pos;
		return At;
	}

	//------------------------------------------1,1--------------------------------------------

std::vector<double> alglib_linear_solve(const arcsim::SpMat<double> &A, const std::vector<double> &b) {
		const int n = b.size();
		alglib::real_2d_array M;
		alglib::real_1d_array x, c;
		M.setlength(n, n);
		x.setlength(n);
		c.setcontent(n, &b[0]);
		for (int i = 0; i < A.m; i++) {
			const arcsim::SpVec<double> &row = A.rows[i];
			for (int j = 0; j < A.n; j++) M(i, j) = 0;
			for (size_t jj = 0; jj < row.indices.size(); jj++) {
				int j = row.indices[jj];
				M(i, j) = row.entries[jj];
			}
		}

		alglib::ae_int_t info;
		alglib::densesolverreport rep;
		alglib::rmatrixsolve(M, n, c, info, rep, x);

		std::vector<double> ret(n);
		for (int i = 0; i < n; i++) ret[i] = x[i];
		return ret;
	}
template <int C>
std::vector<arcsim::Vec<C> > alglib_linear_solve_vec(const arcsim::SpMat<arcsim::Mat<C, C> > &A, const std::vector<arcsim::Vec<C> > &b) {
	const int n = b.size() * C;
	alglib::real_2d_array M;
	alglib::real_1d_array x, c;
	M.setlength(n, n);
	x.setlength(n);
	c.setlength(n);
	for (int i = 0; i < n; i++) c[i] = b[i / C][i % C];
	for (int i = 0; i < A.m; i++) {
		const arcsim::SpVec<arcsim::Mat<C, C> > &row = A.rows[i];
		for (size_t jj = 0; jj < row.indices.size(); jj++) {
			int j = row.indices[jj];
			for (int si = 0; si < C; si++)
				for (int sj = 0; sj < C; sj++) M(i * C + si, j * C + sj) = row.entries[jj](si, sj);
		}
	}

	alglib::ae_int_t info;
	alglib::densesolverreport rep;
	alglib::rmatrixsolve(M, n, c, info, rep, x);

	std::vector<arcsim::Vec<C> > ret(n);
	for (int i = 0; i < n; i++) ret[i / C][i % C] = x[i];

	return ret;
}

	//----------------------------main-------------------------------------------------------------------
	std::vector<double> taucs_linear_solve(const arcsim::SpMat<double> &A, const std::vector<double> &b) {
		if (b.size() < 20) return alglib_linear_solve(A, b);//1-1

															// taucs_logfile("stdout");
		taucs_ccs_matrix *Ataucs = sparse_to_taucs(A);//1-2
		std::vector<double> x(b.size());
		char *options[] = { (char *)"taucs.factor.LLT=true", NULL };
		int retval = taucs_linsolve(Ataucs, NULL, 1, &x[0], (double *)&b[0], options, NULL);
		if (retval != TAUCS_SUCCESS) {
			std::cerr << "Error: TAUCS(scalar) failed with return value " << retval << std::endl;
			//segfault();
			exit(EXIT_FAILURE);
		}
		taucs_ccs_free(Ataucs);
		return x;
	}template <int m>
		std:: vector<arcsim::Vec<m> > taucs_linear_solve(const arcsim::SpMat<arcsim::Mat<m, m> > &A, const std::vector<arcsim::Vec<m> > &b) {
		if (b.size() < 6) return alglib_linear_solve_vec(A, b);

		// taucs_logfile("stdout");
		taucs_ccs_matrix *Ataucs = sparse_to_taucs(A);
		std::vector<arcsim::Vec<m> > x(b.size());
		char *options[] = { (char *)"taucs.factor.LLT=true", NULL };
		int retval = taucs_linsolve(Ataucs, NULL, 1, &x[0], (double *)&b[0], options, NULL);
		if (retval != TAUCS_SUCCESS) {
			std::cerr << "Error: TAUCS(vector) failed with return value " << retval << std::endl;
			arcsim::segfault();
			exit(EXIT_FAILURE);
		}
		taucs_ccs_free(Ataucs);
		return x;
	}

	template std::vector<arcsim::Vec3> taucs_linear_solve(const arcsim::SpMat<arcsim::Mat3x3> &A, const std::vector<arcsim::Vec3> &b);


//}