#pragma once
#include "plasticity.hpp"
#include "mesh_struct.hpp"
#include "geometry.hpp"
#include "sub_struct.hpp"



namespace arcsim{
Mat3x3 stretch_plasticity_from_embedding(const Face *face)
{
	double limit = face->material->plastic_limit;
	Mat3x3 F_ps = derivative(pos<PS>(face->v[0]->node), pos<PS>(face->v[1]->node), pos<PS>(face->v[2]->node), normal<PS>(face), face);
	SVD<3, 3> svd = singular_value_decomposition(F_ps);
	for (int j = 0; j < 3; j++) svd.s[j] = clamp(1.0 / svd.s[j], 1. / limit, limit);
	Mat3x3 M = svd.Vt.t() *diag(svd.s)*svd.Vt;
	return M;
}


Mat3x3 edges_to_face(const Vec3 &theta, const Face *face) {
	Mat3x3 S;
	Vec3 n = normal<MS>(face);
	for (int e = 0; e < 3; e++) {
		// const Edge *edge = face->adje[e];
		Vec3 e_mat = face->v[PREV(e)]->u - face->v[NEXT(e)]->u, t_mat = cross(normalize(e_mat), n);
		S -= 1 / 2. * theta[e] * norm(e_mat) * outer(t_mat, t_mat);
	}
	S /= face->a;
	return S;
}


}