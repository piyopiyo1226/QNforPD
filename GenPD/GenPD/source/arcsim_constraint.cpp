#include <algorithm>
#include "arcsim_constraint.hpp"
#include "magic.hpp"


namespace arcsim {

	double IneqCon::value(int *sign) {
		if (sign) *sign = 1;
		double d = 0;
		for (int i = 0; i < 4; i++)
			if (nodes[i]) d += w[i] * dot(n, nodes[i]->x);
		d -= global_index::magic.repulsion_thickness;
		return d;
	}


	MeshGrad IneqCon::gradient() {
		MeshGrad grad;
		for (int i = 0; i < 4; i++)
			if (nodes[i]) grad.push_back(MeshGradV(nodes[i], w[i] * n));
		return grad;
	}

	MeshGrad IneqCon::project() {
		double d = value() + global_index::magic.repulsion_thickness - global_index::magic.projection_thickness;
		if (d >= 0) return MeshGrad();
		double inv_mass = 0;
		for (int i = 0; i < 4; i++)
			if (nodes[i] && free[i]) inv_mass += sq(w[i]) / nodes[i]->m;
		MeshGrad dx;
		for (int i = 0; i < 4; i++)
			if (nodes[i] && free[i]) dx.push_back(MeshGradV(nodes[i], -(w[i] / nodes[i]->m) / inv_mass * n * d));
		return dx;
	}

	double violation(double value) { return (std::max)(-value, 0.); }
	
	double IneqCon::energy(double value) {
		double v = violation(value);
		return stiff * v * v * v / global_index::magic.repulsion_thickness / 6;
	}
	double IneqCon::energy_grad(double value) { return -stiff * sq(violation(value)) / global_index::magic.repulsion_thickness / 2; }
	double IneqCon::energy_hess(double value) { return stiff * violation(value) / global_index::magic.repulsion_thickness; }
	
	MeshGrad IneqCon::friction(double dt, MeshHess &jac) {
		if (mu == 0) return MeshGrad();
		double fn = abs(energy_grad(value()));
		if (fn == 0) return MeshGrad();
		Vec3 v = Vec3(0);
		double inv_mass = 0;
		for (int i = 0; i < 4; i++) {
			if (nodes[i]) {
				v += w[i] * nodes[i]->v;
				if (free[i]) inv_mass += sq(w[i]) / nodes[i]->m;
			}
		}
		Mat3x3 T = Mat3x3(1) - outer(n, n);
		double vt = norm(T * v);
		double f_by_v = (std::min)(mu * fn / vt, 1 / (dt * inv_mass));
		// double f_by_v = mu*fn/max(vt, 1e-1);
		MeshGrad force;
		for (int i = 0; i < 4; i++) {
			if (nodes[i] && free[i]) {
				force.push_back(MeshGradV(nodes[i], -w[i] * f_by_v * T * v));
				for (int j = 0; j < 4; j++) {
					if (free[j]) {
						jac.push_back(MeshHessV(nodes[i], nodes[j], -w[i] * w[j] * f_by_v * T));
					}
				}
			}
		}
		return force;
	}
	
	bool IneqCon::contains(Node *_node) {
		for (int i = 0; i < 4; i++) {
			if (w[i] && nodes[i] == _node) return true;
		}
		return false;
	}

	double EqCon::value(int *sign) {
		if (sign) *sign = 0;
		return dot(n, node->x - x);
	}
	MeshGrad EqCon::gradient() {
		MeshGrad grad;
		grad.push_back(MeshGradV(node, n));
		return grad;
	}
	MeshGrad EqCon::project() { return MeshGrad(); }
	double EqCon::energy(double value) { return stiff * sq(value) / 2.; }
	double EqCon::energy_grad(double value) { return stiff * value; }
	double EqCon::energy_hess(double value) { return stiff; }
	MeshGrad EqCon::friction(double dt, MeshHess &jac) { return MeshGrad(); }
	bool EqCon::contains(Node *_node) { return node == _node; }




}