
#ifndef NEAROBS_HPP
#define NEAROBS_HPP

#include <map>
#include "mesh_struct.hpp"
#include "collision.hpp"//AccelStruct

#include "geometry.hpp"

namespace arcsim {

	std::map<Node *, Plane> nearest_obstacle_planes(const std::vector<Node *> &nodes,const std::vector<AccelStruct *> &obs_accs);

		struct NearPoint
		{
			double d;
			Vec3 x;
			NearPoint(double d, const Vec3 &v ):d(d),x(x){}
		};
	
		Vec3 nearest_point(const Vec3 &x, const std::vector<AccelStruct*>  &accs, double dmin);
		void update_nearest_point(const Vec3 &x, BVHNode *node, NearPoint &p);

}
#endif