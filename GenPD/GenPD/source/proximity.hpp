#pragma once
#include <vector>
#include "arcsim_constraint.hpp"
#include "magic.hpp"
#include "geometry.hpp"

namespace arcsim {
	template <typename T>
	struct Min {
		double key;
		T val;
		Min()
			: key(infinity)
			, val() {}
		void add(double key, T val) {
			if (key < this->key) {
				this->key = key;
				this->val = val;
			}
		}
	};
	void make_proxy_constraints(Mesh &mesh, CollisionProxy &proxy, std::vector<Constraint *> &cons);
	std::vector<Constraint *> proximity_constraints(std::vector<Mesh *> &meshes, const std::vector<Mesh *> &obs_meshes, double mu, double mu_obs, bool proxy_only, bool obs_only);

	void find_proximities(const Face *face0, const Face *face1);

	void add_proximity(const Node *node, const Face *face);
	void add_proximity(const Edge *edge0, const Edge *edge1);
	void add_proximity(Node *node, Edge *edge) ;
	bool in_wedge(double w, const Edge *edge0, const Edge *edge1);// {//8-1-6

	std::vector<Constraint *> proximity_constraints(std::vector<Mesh *> &meshes, const std::vector<Mesh *> &obs_meshes, double mu,
		double mu_obs, bool proxy_only, bool obs_only);

	Constraint *make_constraint(const Node *node, const Face *face, double mu, double mu_obs);// {//8-2-1
	Constraint *make_constraint(const Edge *edge, const Node *node, double mu, double mu_obs);// {
	Constraint *make_constraint(const Edge *edge0, const Edge *edge1, double mu, double mu_obs);//


	double area_cached(const Face *face);
	double area_cached(const Node *node);//

}