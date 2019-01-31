#ifndef CONSTRAINT_HPP
#define CONSTRAINT_HPP
#include "sub_struct.hpp"

namespace arcsim {

	struct MeshGradV {
		Node *node;
		Vec3 f;
		MeshGradV(Node *i, const Vec3 &f)
			: node(i)
			, f(f) {}
	};
	struct MeshHessV {
		Node *i;
		Node *j;
		Mat3x3 J;
		MeshHessV(Node *i, Node *j, const Mat3x3 &J)
			: i(i)
			, j(j)
			, J(J) {}
	};

	typedef std::vector<MeshGradV> MeshGrad;
	typedef std::vector<MeshHessV> MeshHess;
	struct Constraint {
		virtual ~Constraint() {};
		virtual double value(int *sign = NULL) = 0;
		virtual MeshGrad gradient() = 0;
		virtual MeshGrad project() = 0;
		virtual bool contains(Node *node) = 0;
		// energy function
		virtual double energy(double value) = 0;
		virtual double energy_grad(double value) = 0;
		virtual double energy_hess(double value) = 0;
		// frictional force
		virtual MeshGrad friction(double dt, MeshHess &jac) = 0;

		// useful for debugging
	};



	struct IneqCon : public Constraint//あとでかみチェック
	{ // n . sum(w[i] verts[i]->x) >= 0
		Node *nodes[4];
		double w[4];
		bool free[4];
		Vec3 n;
		double a;//area
		double mu;//friction
		double stiff;
		double value(int *sign = NULL);
		bool contains(Node *node);
		MeshGrad gradient();
		MeshGrad project();
		double energy(double value);
		double energy_grad(double value);//vaue
		double energy_hess(double value);
		MeshGrad friction(double dt, MeshHess &jac);
	
		//useful for debug
	
	//	void serializer(Serialize &s, const std::string &name);
	
	};


	struct EqCon : public Constraint {
		// n . (node->x - x) = 0
		Node *node;
		Vec3 x, n;
		double stiff;
		double value(int *sign = NULL);
		bool contains(Node *node);
		MeshGrad gradient();
		MeshGrad project();
		double energy(double value);
		double energy_grad(double value);
		double energy_hess(double value);
		MeshGrad friction(double dt, MeshHess &jac);
	};

	//------------------------------------------for simulation ---------------------//

}
#endif