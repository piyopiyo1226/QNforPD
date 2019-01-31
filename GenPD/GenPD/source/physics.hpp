
#include "arcsim_constraint.hpp"
#include "Sparce.hpp"
#include <iostream>

namespace arcsim {
	template <Space s>
	double internal_energy(const std::vector<Face *> &faces, const std::vector<Edge *> &edges);

	template <Space s>
	Mat3x3 deformation_gradient(const Face *face);



	template <Space s>
	void add_internal_forces(const std::vector<Face *> &faces, const std::vector<Edge *> &edges, SpMat<Mat3x3> &A,
		std::vector<Vec3> &b, double dt);


	inline Vec3 get_subvec(const double *x, int i){ return Vec3(x[i * 3 + 0], x[i * 3 + 1], x[i * 3 + 2]); }

	double bending_coeff(const Edge *edge, double theta); 
	Mat3x3 material_model(const Face *face, const Mat3x3 &G);
	double constraint_energy(const std::vector<Constraint *> &cons);
	void add_constraint_forces(const std::vector<Constraint *> &cons, SpMat<Mat3x3> &A, std::vector<Vec3> &b, double dt);


}