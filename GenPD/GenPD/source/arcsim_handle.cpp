#include "arcsim_handle.hpp"
#include "magic.hpp"

namespace arcsim
{

	static Vec3 directions[3] = { Vec3(1, 0, 0), Vec3(0, 1, 0), Vec3(0, 0, 1) };
	void add_position_constraints(const Node *node, const Vec3 &x, double stiff, std::vector<Constraint *> &cons);


	Transformation normalize(const Transformation &T) {
		Transformation T1 = T;
		T1.rotation = normalize(T1.rotation);
		return T1;
	}


	std::vector<Constraint *> NodeHandle::get_constraints(double t) {
		double s = strength(t);
		if (!s) return std::vector<Constraint *>();
		if (!activated) {
			//// handle just got started, fill in its original position
			x0 = motion ? inverse(normalize(motion->pos(t))).apply(node->x) : node->x;
			activated = true;
		}
		Vec3 x = motion ? normalize(motion->pos(t)).apply(x0) : x0;
		std::vector<Constraint *> cons;
		add_position_constraints(node, x, s * global_index::magic.handle_stiffness, cons);
		return cons;
	}

	void add_position_constraints(const Node *node, const Vec3 &x, double stiff, std:: vector<Constraint *> &cons) {
		for (int i = 0; i < 3; i++) {
			EqCon *con = new EqCon;
			con->node = (Node *)node;
			con->x = x;
			con->n = directions[i];
			con->stiff = stiff;
			cons.push_back(con);
		}
	}

}