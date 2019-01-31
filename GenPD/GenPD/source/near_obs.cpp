
#include "near_obs.hpp"
#include "magic.hpp"

namespace arcsim {

		
		 std::map<Node *, Plane> nearest_obstacle_planes(const std::vector<Node*> &nodes, const std::vector<AccelStruct * > &obs_accs)
		{
			const double dmin = 10 * global_index::magic.repulsion_thickness;
			std::vector<Plane> planes(nodes.size(), Plane(Vec3(0), Vec3(0)));
			for (int64_t n = 0; n < nodes.size(); n++)
			{
				Vec3 x = nodes[n]->x;
				Vec3 p = nearest_point(x, obs_accs, dmin);//check
				if (p != x) planes[n].x0 = p;
				planes[n].n = normalize(x - p);
			}
			std::map<Node *, Plane> plane_map;
			for (size_t i = 0; i < planes.size(); i++)
				if (norm2(planes[i].n) != 0) plane_map[(Node *)nodes[i]] = planes[i];
			return plane_map;
		}



		Vec3 nearest_point(const Vec3 &x, const std::vector<AccelStruct*>  &accs, double dmin)
		{
			NearPoint p(dmin, x);
			for (int i = 0; i < (int)accs.size(); i++)//addelstructのおおきさは何をいみしてる？
			{
				if (accs[i]->root) update_nearest_point(x, accs[i]->root, p);
				return p.x;
			}

		}



		////-------------------------------for search nearest pint -------------------
		double point_box_dinstance(const Vec3 &x, const BOX &box)
			{
				Vec3 xp = Vec3(clamp(x[0], (double)box._dist[0], (double)box._dist[9]),//TODO: clampドウニュウシテネ
					clamp(x[1], (double)box._dist[1], (double)box._dist[10]),
					 clamp(x[2], (double)box._dist[2], (double)box._dist[11]));
				return norm(x - xp);
		}
		
		
		void update_nearest_point(const Vec3 &x, const Face *face, NearPoint &p)
		{
				
				Vec3 n;
				double w[4];
				double d = unsigned_vf_distance(x, face->v[0]->node->x, face->v[1]->node->x, face->v[2]->node->x, &n, w);
				if (d < p.d)
				{
					p.d = d;
					p.x = -(w[1] * face->v[0]->node->x + w[2] * face->v[1]->node->x + w[3] * face->v[2]->node->x);
				}
		}
		
		void update_nearest_point(const Vec3 &x, BVHNode *node, NearPoint &p)//
		{
				if (node->isLeaf()) update_nearest_point(x, node->getFace(), p);
				else {
					double d = 0;// point_box_dinstance(x, node->_box);
					if (d >= p.d) return;
					update_nearest_point(x, node->getRightChild(), p);
					update_nearest_point(x, node->getLeftChild(), p);
				}
		}


}