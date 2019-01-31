#pragma once
#include <vector>
#include "arcsim_constraint.hpp"
#include "mesh_struct.hpp"

namespace arcsim
{

	struct Handle {
		double start_time, end_time, fade_time;
		virtual ~Handle() {};
		virtual std::vector<Constraint *> get_constraints(double t) = 0;
		virtual std::vector<Node *> get_nodes() = 0;
		bool active(double t) { return t >= start_time && t <= end_time; }
		double strength(double t) {
			if (t < start_time || t > end_time + fade_time) return 0;
			if (t <= end_time) return 1;
			double s = 1 - (t - end_time) / (fade_time + 1e-6);
			return sq(sq(s));
		}
		virtual void add_forces(double t, std::vector<Vec3> &fext, std::vector<Mat3x3> &Jext) {}
	};

	struct NodeHandle : public Handle {
		Node *node;
		const Motion *motion;
		bool activated;
		Vec3 x0;
		NodeHandle()
			: activated(false) {}
		std::vector<Constraint *> get_constraints(double t);
		std::vector<Node *> get_nodes() { return std::vector<Node *>(1, node); }
	};

	void add_position_constraints(const Node *node, const Vec3 &x, double stiff, std::vector<Constraint *> &cons);//



}