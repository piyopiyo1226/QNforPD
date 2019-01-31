#include "physics.hpp"
#include "remesh.hpp"
#include "blockvectors.hpp"

namespace arcsim
{
	namespace global_index{
	static const bool verbose = false;
	}


	// A = dt^2 J + dt damp J
	// b = dt f + dt^2 J v + dt damp J v

	typedef Mat<9, 9> Mat9x9;
	typedef Mat<9, 6> Mat9x6;
	typedef Mat<6, 6> Mat6x6;
	typedef Mat<4, 6> Mat4x6;
	typedef Mat<3, 4> Mat3x4;
	typedef Mat<4, 9> Mat4x9;
	typedef Vec<9> Vec9;

	// A kronecker B = [a11 B, a12 B, ..., a1n B;
	//                  a21 B, a22 B, ..., a2n B;
	//                   ... ,  ... , ...,  ... ;
	//                  am1 B, am2 B, ..., amn B]
	template <int m, int n, int p, int q>
	Mat<m * p, n * q> kronecker(const Mat<m, n> &A, const Mat<p, q> &B) {
		Mat<m * p, n * q> C;
		for (int i = 0; i < m; i++)
			for (int j = 0; j < n; j++)
				for (int k = 0; k < p; k++)
					for (int l = 0; l < q; l++) C(i * p + k, j * q + l) = A(i, j) * B(k, l);
		return C;
	}

	template <int n>
	Mat<1, n> rowmat(const Vec<n> &v) {
		Mat<1, n> A;
		for (int i = 0; i < n; i++) A(0, i) = v[i];
		return A;
	}



	template <Space s>
	double bending_energy(const Edge *edge) {
		if (!edge->adjf[0] || !edge->adjf[1]) return 0;

		double theta = dihedral_angle<s>(edge);
		return bending_coeff(edge, theta) * sq(theta - edge->theta_ideal) / 4;
	}
	template <Space s>
	double stretching_energy(const Face *face) {
		Mat3x3 F = deformation_gradient<WS>(face);
		Mat3x3 G = (F.t() * F - Mat3x3(1)) * 0.5;
		Mat3x3 sigma = material_model(face, G);

		return face->a * 0.5 * inner(sigma, G);
	}


	template <Space s>
	double internal_energy(const std::vector<Face *> &faces, const std::vector<Edge *> &edges) {
		double E = 0;
		for (size_t f = 0; f < faces.size(); f++) E += stretching_energy<s>(faces[f]);
		for (size_t e = 0; e < edges.size(); e++) {
			E += bending_energy<s>(edges[e]);
		}
		return E;
	}
	template double internal_energy<PS>(const std::vector<Face *> &, const std::vector<Edge *> &);
	template double internal_energy<WS>(const std::vector<Face *> &, const std::vector<Edge *> &);

	Vec<3, int> indices(const Node *n0, const Node *n1, const Node *n2) {
		Vec<3, int> ix;
		ix[0] = n0->active() ? n0->index : -1;
		ix[1] = n1->active() ? n1->index : -1;
		ix[2] = n2->active() ? n2->index : -1;
		return ix;
	}

	Vec<4, int> indices(const Node *n0, const Node *n1, const Node *n2, const Node *n3) {
		Vec<4, int> ix;
		ix[0] = n0->active() ? n0->index : -1;
		ix[1] = n1->active() ? n1->index : -1;
		ix[2] = n2->active() ? n2->index : -1;
		ix[3] = n3->active() ? n3->index : -1;
		return ix;
	}


	template <Space s>
	Mat3x3 deformation_gradient(const Face *face) {//21-1
		return derivative(pos<s>(face->v[0]->node), pos<s>(face->v[1]->node), pos<s>(face->v[2]->node), normal<s>(face),
			face) *
			face->Sp_str;
	}


	template <Space s>
	void add_internal_forces(const std::vector<Face *> &faces, const std::vector<Edge *> &edges, SpMat<Mat3x3> &A, std::vector<Vec3> &b,
		double dt) {
		for (size_t f = 0; f < faces.size(); f++) {
			const Face *face = faces[f];
			const Node *n0 = face->v[0]->node, *n1 = face->v[1]->node, *n2 = face->v[2]->node;
			Vec9 vs = mat_to_vec(Mat3x3(n0->v, n1->v, n2->v));
			std::pair<Mat9x9, Vec9> membF = stretching_force<s>(face);
			Mat9x9 J = membF.first;
			Vec9 F = membF.second;
			if (dt == 0) {
				add_submat(-J, indices(n0, n1, n2), A);
				add_subvec(F, indices(n0, n1, n2), b);
			}
			else {
				double damping = face->material->damping;
				add_submat(-dt * (dt + damping) * J, indices(n0, n1, n2), A);
				add_subvec(dt * (F + (dt + damping) * J * vs), indices(n0, n1, n2), b);
			}
		}
		for (size_t e = 0; e < edges.size(); e++) {
			const Edge *edge = edges[e];
			if (!edge->adjf[0] || !edge->adjf[1]) continue;
			std::pair<Mat12x12, Vec12> bendF = bending_force<s>(edge);
			const Node *n0 = edge->n[0], *n1 = edge->n[1], *n2 = edge_opp_vert(edge, 0)->node,
				*n3 = edge_opp_vert(edge, 1)->node;
			Vec12 vs = mat_to_vec(Mat3x4(n0->v, n1->v, n2->v, n3->v));
			Mat12x12 J = bendF.first;
			Vec12 F = bendF.second;

			if (dt == 0) {
				add_submat(-J, indices(n0, n1, n2, n3), A);
				add_subvec(F, indices(n0, n1, n2, n3), b);
			}
			else {
				double damping = (edge->adjf[0]->material->damping + edge->adjf[1]->material->damping) * 0.5;
				add_submat(-dt * (dt + damping) * J, indices(n0, n1, n2, n3), A);
				add_subvec(dt * (F + (dt + damping) * J * vs), indices(n0, n1, n2, n3), b);
			}
		}
	}
	template void add_internal_forces<PS>(const std::vector<Face *> &, const std::vector<Edge *> &, SpMat<Mat3x3> &, std::vector<Vec3> &,
		double);
	template void add_internal_forces<WS>(const std::vector<Face *> &, const std::vector<Edge *> &, SpMat<Mat3x3> &, std::vector<Vec3> &,
		double);



	template <Space s>
	std::pair<Mat9x9, Vec9> stretching_force(const Face *face) {
		const Material *mat = face->material;
		// compute stress, strain
		const Vec3 x[3] = { pos<s>(face->v[0]->node), pos<s>(face->v[1]->node), pos<s>(face->v[2]->node) };
		Mat3x3 F = deformation_gradient<s>(face);
		Mat3x3 G = (F.t() * F - Mat3x3(1)) * 0.5;
		double weakening_mult = 1 / (1 + mat->weakening * face->damage);

		Mat3x3 Y = face->invDm * face->Sp_str;
		Mat3x3 D = Mat3x3::rows(-Y.row(0) - Y.row(1), Y.row(0), Y.row(1));
		Mat<3, 9> DD[3] = { kronecker(rowmat(D.col(0)), Mat3x3(1)), kronecker(rowmat(D.col(1)), Mat3x3(1)),
			kronecker(rowmat(D.col(2)), Mat3x3(1)) };
		Vec9 X;
		for (int i = 0; i < 9; i++) X[i] = x[i / 3][i % 3];
		Vec3 f[3] = { DD[0] * X, DD[1] * X, DD[2] * X };

		Vec9 grad_f(0);
		Mat9x9 hess_f(0);
		if (mat->use_dde) {
			Vec4 k = stretching_stiffness(reduce_xy(G), mat->dde_stretching) * weakening_mult;

			const Mat<3, 9> &Du = DD[0];
			const Mat<3, 9> &Dv = DD[1];
			Vec9 fuu = Du.t() * f[0], fvv = Dv.t() * f[1], fuv = (Du.t() * f[1] + Dv.t() * f[0]) / 2.;
			grad_f = k[0] * G(0, 0) * fuu + k[2] * G(1, 1) * fvv + k[1] * (G(0, 0) * fvv + G(1, 1) * fuu) +
				2 * k[3] * G(0, 1) * fuv;
			hess_f = k[0] * (outer(fuu, fuu) + (std::max)(G(0, 0), 0.) * Du.t() * Du) +
				k[2] * (outer(fvv, fvv) + (std::max)(G(1, 1), 0.) * Dv.t() * Dv) +
				k[1] * (outer(fuu, fvv) + (std::max)(G(0, 0), 0.) * Dv.t() * Dv + outer(fvv, fuu) +
				(std::max)(G(1, 1), 0.) * Du.t() * Du) +
				2. * k[3] * (outer(fuv, fuv));
			// ignoring G(0,1)*(Du.t()*Dv+Dv.t()*Du)/2. term
			// because may not be positive definite
			return std::make_pair(-face->a * hess_f, -face->a * grad_f);
		}
		else {
			double gf = -face->a * mat->alt_stretching;
			// Vec9 d_trace = 0.5 * (DD[0].t()*DD[0]+DD[1].t()*DD[1]+DD[2].t()*DD[2]) * X;
			Mat3x3 Gc = max(G, Mat3x3(0));  // posdef
			for (int i = 0; i < 3; i++)
				for (int j = 0; j <= i; j++) {
					Vec9 dG = 0.5 * (DD[i].t() * f[j] + DD[j].t() * f[i]);
					Mat9x9 dG2 = 0.5 * (DD[i].t() * DD[j] + DD[j].t() * DD[i]);
					Vec9 d_trace_c = 0.5 * dG;  // posdef

					if (i == j)
						grad_f += gf * (1.0 - mat->alt_poisson) * G(i, j) * dG + gf * mat->alt_poisson * trace(G) * dG;
					else
						grad_f += 2.0 * gf * (1.0 - mat->alt_poisson) * G(i, j) * dG;

					if (i == j)
						hess_f += gf * (1.0 - mat->alt_poisson) * (outer(dG, dG) + Gc(i, j) * dG2) +
						gf * mat->alt_poisson * (trace(Gc) * dG2 + outer(d_trace_c, dG));
					else
						hess_f += 2.0 * gf * (1.0 - mat->alt_poisson) * (outer(dG, dG) + Mat9x9(0));  // posdef
				}
			return std::make_pair(hess_f, grad_f);
		}
	}

	template <int m, int n>
	Mat<3, 3> submat3(const Mat<m, n> &A, int i, int j) {
		Mat3x3 Asub;
		for (int k = 0; k < 3; k++)
			for (int l = 0; l < 3; l++) Asub(k, l) = A(i * 3 + k, j * 3 + l);
		return Asub;
	}
	template <int n>
	Vec<3> subvec3(const Vec<n> &b, int i) {
		Vec3 bsub;
		for (int k = 0; k < 3; k++) bsub[k] = b[i * 3 + k];
		return bsub;
	}

	template <int m>
	void add_subvec(const Vec<m * 3> &bsub, const Vec<m, int> &ix, std::vector<Vec3> &b) {
		for (int i = 0; i < m; i++) {
			if (ix[i] < 0) continue;
			b[ix[i]] += subvec3(bsub, i);
		}
	}

	template <int m>
	void add_submat(const Mat<m * 3, m * 3> &Asub, const Vec<m, int> &ix, SpMat<Mat3x3> &A) {
		for (int i = 0; i < m; i++) {
			if (ix[i] < 0) continue;
			for (int j = 0; j < m; j++) {
				if (ix[j] < 0) continue;
				A(ix[i], ix[j]) += submat3(Asub, i, j);
			}
		}
	}

	double distance(const Vec3 &x, const Vec3 &a, const Vec3 &b) {
		Vec3 e = b - a;
		Vec3 xp = e * dot(e, x - a) / dot(e, e);
		// return norm((x-a)-xp);
		return (std::max)(norm((x - a) - xp), 1e-3 * norm(e));
	} 
	typedef Mat<12, 12> Mat12x12;
	typedef Vec<12> Vec12;

	double bending_coeff(const Edge *edge, double theta) {
		const Face *face0 = edge->adjf[0], *face1 = edge->adjf[1];
		double a = face0->a + face1->a;
		double l = norm(edge->n[1]->x - edge->n[0]->x);
		double ke0 = (face0->material->use_dde) ? bending_stiffness(edge, 0, face0->material->dde_bending, l, theta)
			: face0->material->alt_bending;
		double ke1 = (face1->material->use_dde) ? bending_stiffness(edge, 1, face1->material->dde_bending, l, theta)
			: face1->material->alt_bending;

		double ke = (std::min)(ke0, ke1);
		double weakening = (std::max)(face0->material->weakening, face1->material->weakening);
		ke *= 1 / (1 + weakening * edge->damage);
		double shape = sq(l) / (2 * a);

		return ke * shape;
	}

	Vec2 barycentric_weights(const Vec3 &x, const Vec3 &a, const Vec3 &b) {
		Vec3 e = b - a;
		double t = dot(e, x - a) / dot(e, e);
		return Vec2(1 - t, t);
	}

	template <Space s>
	std::pair<Mat12x12, Vec12> bending_force(const Edge *edge) {
		const Face *face0 = edge->adjf[0], *face1 = edge->adjf[1];
		if (!face0 || !face1) return std::make_pair(Mat12x12(0), Vec12(0));
		double theta = dihedral_angle<s>(edge);
		Vec3 x0 = pos<s>(edge->n[0]), x1 = pos<s>(edge->n[1]), x2 = pos<s>(edge_opp_vert(edge, 0)->node),
			x3 = pos<s>(edge_opp_vert(edge, 1)->node);
		double h0 = distance(x2, x0, x1), h1 = distance(x3, x0, x1);
		Vec3 n0 = normal<s>(face0), n1 = normal<s>(face1);
		Vec2 w_f0 = barycentric_weights(x2, x0, x1), w_f1 = barycentric_weights(x3, x0, x1);
		Vec12 dtheta = mat_to_vec(
			Mat3x4(-(w_f0[0] * n0 / h0 + w_f1[0] * n1 / h1), -(w_f0[1] * n0 / h0 + w_f1[1] * n1 / h1), n0 / h0, n1 / h1));

		double coeff = bending_coeff(edge, theta);
		return std::make_pair(-coeff * outer(dtheta, dtheta) / 2., -coeff * (theta - edge->theta_ideal) * dtheta / 2.);
	}


	Mat3x3 material_model(const Face *face, const Mat3x3 &G) {
		const Material *mat = face->material;
		double weakening_mult = 1 / (1 + mat->weakening * face->damage);
		if (mat->use_dde) {
			Vec4 k = stretching_stiffness(reduce_xy(G), mat->dde_stretching) * weakening_mult;
			Mat3x3 sigma(Vec3(k[0] * G(0, 0) + k[1] * G(1, 1), 0.5 * k[3] * G(0, 1), 0),
				Vec3(0.5 * k[3] * G(1, 0), k[2] * G(1, 1) + k[1] * G(0, 0), 0), Vec3(0, 0, 0));
			return sigma;
		}
		else {
			double A = mat->alt_stretching * weakening_mult;
			Mat3x3 sigma = A * (1.0 - mat->alt_poisson) * G + Mat3x3(A * mat->alt_poisson * trace(G));
			return sigma;
		}
	}
	void add_constraint_forces(const std::vector<Constraint *> &cons, SpMat<Mat3x3> &A, std::vector<Vec3> &b, double dt) {
		for (int c = 0; c < (int)cons.size(); c++) {
			double value = cons[c]->value();
			double g = cons[c]->energy_grad(value);
			double h = cons[c]->energy_hess(value);
			MeshGrad grad = cons[c]->gradient();
			// f = -g*grad
			// J = -h*outer(grad,grad)
			double v_dot_grad = 0;
			for (MeshGrad::iterator it = grad.begin(); it != grad.end(); it++) {
				v_dot_grad += dot(it->f, it->node->v);
			}
			for (MeshGrad::iterator it = grad.begin(); it != grad.end(); it++) {
				const Node *nodei = it->node;
				if (!nodei->active()) continue;
				int ni = nodei->index;
				if (ni < 0 || ni >= A.m) continue;

				for (MeshGrad::iterator jt = grad.begin(); jt != grad.end(); jt++) {
					const Node *nodej = jt->node;
					if (!nodej->active()) continue;
					int nj = nodej->index;
					if (nj < 0 || nj >= A.n) continue;

					if (dt == 0)
						A(ni, nj) += h * outer(it->f, jt->f);
					else
						A(ni, nj) += dt * dt * h * outer(it->f, jt->f);
				}
				if (dt == 0)
					b[ni] -= g * it->f;
				else
					b[ni] -= dt * (g + dt * h * v_dot_grad) * it->f;
			}
		}
	}

	double constraint_energy(const std::vector<Constraint *> &cons) {
		double E = 0;
		for (int c = 0; c < (int)cons.size(); c++) {
			double value = cons[c]->value();
			double e = cons[c]->energy(value);
			E += e;
		}
		return E;
	}


}