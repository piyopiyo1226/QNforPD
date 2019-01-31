#include "remesh.hpp"
#include "frac.hpp"
#include "optim.hpp"
#include "arcsim_simulation.hpp"
#include "magic.hpp"
#include "tensor.hpp"
#include "referenceshape.hpp"
#include "proximity.hpp"
#include "plasticity.hpp"
#include "physics.hpp"
#include <algorithm>
#include <math.h>

namespace arcsim {

		Arc_sim sim;

		namespace global_index {
			static const bool verbose = false;
		}



		RemeshOp RemeshOp::inverse() const {
			RemeshOp iop;
			iop.added_verts = removed_verts;
			iop.removed_verts = added_verts;
			iop.added_nodes = removed_nodes;
			iop.removed_nodes = added_nodes;
			iop.added_edges = removed_edges;
			iop.removed_edges = added_edges;
			iop.added_faces = removed_faces;
			iop.removed_faces = added_faces;
			return iop;
		}


		void RemeshOp::cancel() {
			delete_all(added_verts);
			delete_all(added_nodes);
			delete_all(added_edges);
			delete_all(added_faces);
			added_edges.clear();
			added_faces.clear();
			added_nodes.clear();
			added_verts.clear();
			removed_edges.clear();
			removed_faces.clear();
			removed_nodes.clear();
			removed_verts.clear();
		} 

		void RemeshOp::apply(Mesh &mesh) const {
			remove_all(removed_faces, mesh);
			remove_all(removed_edges, mesh);
			remove_all(removed_nodes, mesh);
			remove_all(removed_verts, mesh);
			add_all(added_verts, mesh);
			add_all(added_nodes, mesh);
			add_all(added_edges, mesh);
			add_all(added_faces, mesh);
			for (size_t f = 0; f < added_faces.size(); f++) compute_ms_data(added_faces[f]);
			for (size_t f = 0; f < added_faces.size(); f++)
				for (int i = 0; i < 3; i++) compute_ms_data(added_faces[f]->v[i]->node);
		}

	void RemeshOp::done() const {
			delete_all(removed_verts);
			delete_all(removed_nodes);
			delete_all(removed_edges);
	}

	void RemeshOp::update(std::vector<Face *> &v) {
		exclude_all(removed_faces, v);
		include_all(added_faces, v);
	}

	void RemeshOp::update(std::vector<Edge *> &v) {
		exclude_all(removed_edges, v);
		include_all(added_edges, v);
	}

	void RemeshOp::update(std::vector<Node *> &v) {
		exclude_all(removed_nodes, v);
		include_all(added_nodes, v);
	}

		void RemeshOp::set_null(std::vector<Edge *> &v) {
			for (size_t i = 0; i < v.size(); i++)
				if (is_in(v[i], removed_edges)) v[i] = 0;
		}

		struct PlasticityStash {//6-2
			std::vector<Mat3x3> Sp_str;
			std::vector<double> theta_ideal;
			PlasticityStash(std::vector<Face *> &faces, std::vector<Edge *> &edges)
				: Sp_str(faces.size())
				, theta_ideal(edges.size()) {
				for (size_t f = 0; f < faces.size(); f++) {
					Face *face = faces[f];
					Sp_str[f] = face->Sp_str;
					face->Sp_str = Mat3x3(1);
				}
				for (size_t e = 0; e < edges.size(); e++) {
					Edge *edge = edges[e];
					theta_ideal[e] = edge->theta_ideal;
					edge->theta_ideal = dihedral_angle<MS>(edge);
				}
			}
			void apply(std::vector<Face *> &faces, std::vector<Edge *> &edges) {
				for (size_t f = 0; f < faces.size(); f++) {
					faces[f]->Sp_str = Sp_str[f];
				}
				for (size_t e = 0; e < edges.size(); e++) {
					edges[e]->theta_ideal = theta_ideal[e];
				}
			}
		};

		void dynamic_remesh(Mesh &mesh, const std::map<Node *, Plane> &planes) {
			delete_spaced_out(mesh);//1
			create_vert_sizing(mesh.verts, planes);//2
			std::vector<Face *> active_faces = mesh.faces;
			//// cout << "before\n"; wait_key();
			flip_edges(0, active_faces, 0, 0);//3
			while (split_worst_edge(0, mesh.edges));//4
			//	;
			//// cout << "post split\n"; wait_key();
			//active_faces = mesh.faces;
			while (improve_some_face(0, active_faces));//10
				//	;
				//// cout << "post collapse\n"; wait_key();
				//compute_ms_data(mesh);//11
		}



	//-----------------------------------1--------------------------------------------------
	void delete_spaced_out(Mesh &mesh) {
			RemeshOp op;
			for (size_t i = 0; i < mesh.verts.size(); i++)
				if (norm2(mesh.verts[i]->node->x) > 1e6) op.removed_verts.push_back(mesh.verts[i]);
	
			if (op.removed_verts.empty()) return;//ここはなんじゃらほい
	
			bool grow = true;
			//while (grow) {
			//	grow = false;
			//	for (size_t i = 0; i < op.removed_verts.size(); i++) {
			//		Vert *v = op.removed_verts[i];
			//		for (size_t j = 0; j < v->adjf.size(); j++) {
			//			Face *f = v->adjf[j];
			//			if (find(f, op.removed_faces) < 0) {
			//				grow = true;
			//				op.removed_faces.push_back(f);
			//				for (int k = 0; k < 3; k++) {
			//					include(f->adje[k], op.removed_edges);//include is in tuil.hpp入っているかどうか判定してはいっていなかったらさいこうびにいれるようにする
			//					include(f->v[k], op.removed_verts);
			//					include(f->v[k]->node, op.removed_nodes);
			//				}
			//			}
			//		}
			//	}
			//}
			std::cout << "deleted " << op.removed_nodes.size() << " spaced out nodes" << std::endl;
			op.apply(mesh);
			op.done();
		}
		//--------------2  main--//
	
		void  create_vert_sizing(std::vector<Vert *> &verts, const std::map<Node *, Plane> &planes) {//
			Remeshing &remeshing = verts[0]->node->mesh->parent->remeshing;
			//ここまチェックした
			std::map<Face *, Mat3x3> face_sizing;
			for (size_t i = 0; i < verts.size(); i++) {
				Mat3x3 sizing(0);
				Vert *vert = verts[i];
				double wsum = 0;
				for (size_t f = 0; f < vert->adjf.size(); f++) {
					Face *face = vert->adjf[f];
					if (!face) continue;
					if (face_sizing.find(face) == face_sizing.end())
						face_sizing[face] = compute_face_sizing(remeshing, face, planes);// 2-1
					sizing += face->a * face_sizing[face];
					wsum += face->a;
				}
				vert->sizing = sizing / wsum;
			}
		}
	
	//----------------------------------end 2-----------------------------------------------
	//--------------------------------------3---------------------------------//
		

		Mat3x3 compute_face_sizing(Remeshing &remeshing, const Face *face, const std::map<Node *, Plane> &planes, bool debug) {
			// project to in-plane 2D
			Mat3x3 base = local_base(normal<MS>(face));
			Mat3x2 UV(base.col(0), base.col(1));
			Mat2x3 UVt = UV.t();

			// Mat2x2 Sp = projected_curvature<PS>(face,UVt);
			Mat2x2 Sw1 = projected_curvature<WS>(face, UVt);
			Mat3x3 Sw2 = derivative(face->v[0]->node->n, face->v[1]->node->n, face->v[2]->node->n, Vec3(0), face);
			// Mat2x2 Mcurvp = !sim.enabled[Simulation::Plasticity] ? Mat2x2(0)
			//              : (Sp.t()*Sp)/sq(remeshing.refine_angle);
			Mat2x2 Mcurvw1 = (Sw1.t() * Sw1) / sq(remeshing.refine_angle);
			Mat2x2 Mcurvw2 = UVt * (Sw2.t() * Sw2) * UV / sq(remeshing.refine_angle);
			Mat3x3 V = derivative(face->v[0]->node->v, face->v[1]->node->v, face->v[2]->node->v, Vec3(0), face);
			Mat2x2 Mvel = UVt * (V.t() * V) * UV / sq(remeshing.refine_velocity);

			Mat2x2 Mcomp = compression_metric(face, Sw2.t() * Sw2, UV, remeshing.refine_compression);
			Mat2x2 Mobs = (planes.empty()) ? Mat2x2(0) : UVt * obstacle_metric(face, planes) * UV;

			Mat2x2 Mfrac = fracture_metric(remeshing, face) / sq(remeshing.size_min);//hoshi

			std::vector<Mat2x2> Ms(6);
			// Ms[0] = Mcurvp;
			Ms[0] = Mcurvw1;
			Ms[1] = Mcurvw2;
			Ms[2] = Mvel;
			Ms[3] = Mcomp;
			Ms[4] = Mobs;
			Ms[5] = Mfrac;

			Mat2x2 s = global_index::magic.combine_tensors ? tensor_max(Ms) : Ms[0] + Ms[1] + Ms[2] + Ms[3] + Ms[4] + Ms[5] + Ms[6];
			// specific resolution request ?
			for (int i = 0; i < 3; i++) {
				if (face->v[0]->node->flag & Node::FlagResolveUni) s = Mat2x2(1.f / sq(remeshing.size_uniform));
				if (face->v[0]->node->flag & Node::FlagResolveMax) s = Mat2x2(1.f / sq(remeshing.size_min));
			}

			if (debug) {
				std::cout << "curvw1 " << norm_F(Mcurvw1) << std::endl;
				std::cout << "curvw2 " << norm_F(Mcurvw2) << std::endl;
				std::cout << "vel " << norm_F(Mvel) << std::endl;
				std::cout << "comp " << norm_F(Mcomp) << std::endl;
				std::cout << "obs " << norm_F(Mobs) << std::endl;
				std::cout << "frac " << norm_F(Mfrac) << std::endl;
				std::cout << "total " << norm_F(s) << std::endl;
			}

			Eig<2> eig = eigen_decomposition(s);
			for (int i = 0; i < 2; i++) {
				eig.l[i] = clamp(eig.l[i], 1.f / sq(remeshing.size_max), 1.f / sq(remeshing.size_min));
			}
			double lmax = (std::max)(eig.l[0], eig.l[1]);
			double lmin = lmax * sq(remeshing.aspect_min);
			for (int i = 0; i < 2; i++)
				if (eig.l[i] < lmin) eig.l[i] = lmin;
			s = eig.Q * diag(eig.l) * eig.Q.t();
			return UV * s * UVt;  // reproject to 3D
		}

		template <int n>
		Mat<n, n> sqrt(const Mat<n, n> &A) {
			Eig<n> eig = eigen_decomposition(A);
			for (int i = 0; i < n; i++) eig.l[i] = eig.l[i] >= 0 ? sqrt(eig.l[i]) : -sqrt(-eig.l[i]);
			return eig.Q * diag(eig.l) * eig.Q.t();
		}

	

		Mat2x2 perp(const Mat2x2 &A) { return Mat2x2(Vec2(A(1, 1), -A(1, 0)), Vec2(-A(0, 1), A(0, 0))); }


		Mat2x2 compression_metric(const Face *face, const Mat3x3 &S2, const Mat3x2 &UV, double c) {
			Mat3x3 F = deformation_gradient<WS>(face);//21-1
			Mat3x3 G = F.t() * F - Mat3x3(1);
			Mat2x2 e = UV.t() * G * UV;
			Mat2x2 e2 = UV.t() * G.t() * G * UV;
			Mat2x2 Sw2 = UV.t() * S2 * UV;

			Mat2x2 D = e2 - 4.0 * sq(c) * perp(Sw2) * global_index::magic.rib_stiffening;//invecotrs
			return get_positive(-e + sqrt(D)) / (2.0 * sq(c));
		}

		Mat3x3 obstacle_metric(const Face *face, const std::map<Node *, Plane> &planes) {
			Mat3x3 o(0);
			for (int v = 0; v < 3; v++) {
				std::map<Node *, Plane>::const_iterator it = planes.find(face->v[v]->node);
				if (it == planes.end()) continue;
				double h[3];
				for (int v1 = 0; v1 < 3; v1++) h[v1] = dot(face->v[v1]->node->x - it->second.x0, it->second.n);
				Vec3 dh = derivative(h[0], h[1], h[2], 0, face);

				o += outer(dh, dh) / sq(h[v]);
			}
			return o / 3.;
		}


		Mat2x2 fracture_metric(Remeshing &remeshing, const Face *face) {//hoshi
			if (remeshing.refine_fracture == 0 || !sim.enabled[Arc_sim::Fracture]) return Mat2x2(0);
			double fmax = 0;
			for (int i = 0; i < 3; i++) {
				Face *f = adj_face(face, i);
				if (f) {
					Mat3x3 sig = compute_sigma(face);
					Vec3 l = eigen_values(sig);
					fmax = (std::max)(l[0], fmax);
				}
			}

			double refine = remeshing.refine_fracture * 0.5 * face->material->toughness;
			double ramp0 = refine * 0.5, ramp1 = refine / 0.5;
			double v = clamp((fmax - ramp0) / (ramp1 - ramp0), 0.0, 1.0);
			return Mat2x2(v);
		}
		//--------------------------flip----------------------//
			bool independent(const Edge *edge, const std::vector<Edge *> &edges) {//3-3
				for (int i = 0; i < (int)edges.size(); i++) {
					Edge *edge1 = edges[i];
					if (edge->n[0] == edge1->n[0] || edge->n[0] == edge1->n[1] || edge->n[1] == edge1->n[0] ||
						edge->n[1] == edge1->n[1])
						return false;
				}
				return true;
			}


			std::vector<Edge *> independent_edges(const std::vector<Edge *> &edges) { //3-2
				std::vector<Edge *> iedges;
				for (int e = 0; e < (int)edges.size(); e++)
					if (independent(edges[e], iedges)) iedges.push_back(edges[e]);//3-3
				return iedges;
			}

			RemeshOp flip_edge(Edge *edge) {
				RemeshOp op;
				Vert *vert0 = edge_vert(edge, 0, 0), *vert1 = edge_vert(edge, 1, 1), *vert2 = edge_opp_vert(edge, 0),
					*vert3 = edge_opp_vert(edge, 1);
				Face *face0 = edge->adjf[0], *face1 = edge->adjf[1];
				double A = area(face0), B = area(face1);
				Mat3x3 sp = (A * face0->Sp_str + B * face1->Sp_str) / (A + B);
				Mat3x3 sb = (A * face0->Sp_bend + B * face1->Sp_bend) / (A + B);
				double damage = (A * face0->damage + B * face1->damage) / (A + B);
				op.removed_edges.push_back(edge);
				op.added_edges.push_back(new Edge(vert2->node, vert3->node, -edge->theta_ideal, edge->preserve));
				op.removed_faces.push_back(face0);
				op.removed_faces.push_back(face1);
				op.added_faces.push_back(new Face(vert0, vert3, vert2, sp, sb, face0->material, damage));
				op.added_faces.push_back(new Face(vert1, vert2, vert3, sp, sb, face1->material, damage));
				embedding_from_plasticity(op.removed_faces);
				op.apply(*edge->n[0]->mesh);
				plasticity_from_embedding(op.added_faces);
				local_pop_filter(op.added_faces);
				return op;
			}


		bool flip_some_edges(MeshSubset *subset, std::vector<Face *> &active_faces,std::vector<Edge *> *update_edges,
				std::vector<Face *> *update_faces) {//flip_some_edge___
				static int n_edges_prev = 0;
				std::vector<Edge *> edges = independent_edges(find_edges_to_flip(active_faces));//3-2, 3-4
				if ((int)edges.size() == n_edges_prev)  // probably infinite loop
					return false;
		
				bool did_flip = false;
				n_edges_prev = edges.size();
				for (size_t e = 0; e < edges.size(); e++) {
					RemeshOp op = flip_edge(edges[e]);//flip_some_edge___1
					if (op.empty()) continue;
		
					did_flip = true;
					if (subset) op.update(subset->active_nodes);
					if (update_edges) op.set_null(*update_edges);
					if (update_faces) op.update(*update_faces);
					op.update(active_faces);
					op.done();
				}
				return did_flip;
			}
		

		
		void flip_edges(MeshSubset *subset, std::vector<Face *> &active_faces, std::vector<Edge *> *update_edges,
				std::vector<Face *> *update_faces) {
				int N = 3 * active_faces.size();
				for (int i = 0; i < N; i++) {  // don't loop without bound
					if (!flip_some_edges(subset, active_faces, update_edges, update_faces)) return;
				}
			}

			std::vector<Edge *> find_edges_to_flip(std::vector<Face *> &active_faces) {//3-4
				std::vector<Edge *> edges;
				for (size_t i = 0; i < active_faces.size(); i++)
					for (int j = 0; j < 3; j++) include(active_faces[i]->adje[j], edges);//utils
				std::vector<Edge *> fedges;
				for (size_t e = 0; e < edges.size(); e++) {
					Edge *edge = edges[e];
					if (is_seam_or_boundary(edge) || edge->preserve != 0 || !should_flip(edge)) continue;//in util //should flip 3-5
					fedges.push_back(edge);
				}
				return fedges;
			}
			
	// from Bossen and Heckbert 1996
				inline bool should_flip2(const Vec2 &x, const Vec2 &y, const Vec2 &z, const Vec2 &w, const Mat2x2 &M) {//3-6
					double area0 = fabs(wedge(z - y, x - y)), area1 = fabs(wedge(x - w, z - w));
					return area0 * dot(x - w, M * (z - w)) + dot(z - y, M * (x - y)) * area1 <
						-global_index::magic.edge_flip_threshold * (area0 + area1);
				}


				bool should_flip(const Edge *edge) {//3-5
					const double max_angle = 40.0 * M_PI / 180.0;
			
					const Vert *vert0 = edge_vert(edge, 0, 0), *vert1 = edge_vert(edge, 0, 1), *vert2 = edge_opp_vert(edge, 0),
						*vert3 = edge_opp_vert(edge, 1);//inline
					Vec3 x = vert0->u, z = vert1->u, w = vert2->u, y = vert3->u;
			
					// don't flip if high angles are involved
					if (fabs(dihedral_angle<WS>(edge)) >= max_angle) return false;//in geometry check 11/20
			
					const Vec3 &n0 = normalize(cross(y - x, w - x)), n1 = normalize(cross(w - z, y - z));
					if (fabs(dihedral_angle(w, y, n0, n1)) >= max_angle) return false;
			
					// don't flip if mean normal is inverted
					if (dot(n0 + n1, normal<MS>(edge->adjf[0]) + normal<MS>(edge->adjf[1])) <= 0) return false;
			
					// project onto a 2D plane which conserves diagonal lengths
					Vec3 u = normalize(x - z), v = normalize((w - y) - dot(w - y, u) * u);
					Mat3x2 A(u, v);
					Mat2x3 At = A.t();
			
					Mat3x3 M = (vert0->sizing + vert1->sizing + vert2->sizing + vert3->sizing) / 4.0;
			
					Mat2x2 Mr = At * M * A;
					return should_flip2(At * x, At * y, At * z, At * w, Mr);///3-6
				}

				//----------------------------------------split-----------------

				double edge_metric (const Vert *vert0, const Vert *vert1) { //4-3
						if (!vert0 || !vert1) return 0;
						Vec3 du = vert0->u - vert1->u;
						return sqrt((dot(du, vert0->sizing * du) + dot(du, vert1->sizing * du)) / 2.);
				}


				double edge_metric(const Edge *edge) {//4-2
						return (std::max)(edge_metric(edge_vert(edge, 0, 0), edge_vert(edge, 0, 1)),
							edge_metric(edge_vert(edge, 1, 0), edge_vert(edge, 1, 1)));//4-3. 4-4(in inlinemesh)
				}



				bool split_worst_edge(MeshSubset *subset, const std::vector<Edge *> &edges) {
						std::vector<Edge *> bad_edges = find_bad_edges(edges);//4-1
						for (size_t e = 0; e < bad_edges.size(); e++) {
							Edge *edge = bad_edges[e];
							if (!edge) continue;
							Node *node0 = edge->n[0], *node1 = edge->n[1];
							RemeshOp op = split_edge(edge, 0.5);//4-5
							for (size_t v = 0; v < op.added_verts.size(); v++) {
								Vert *vertnew = op.added_verts[v];
								Vert *v0 = adjacent_vert(node0, vertnew), *v1 = adjacent_vert(node1, vertnew);
								vertnew->sizing = 0.5 * (v0->sizing + v1->sizing);
							}
							if (subset) op.update(subset->active_nodes);
							op.set_null(bad_edges);
							op.done();
							if (global_index::verbose) std::cout << "Split " << node0 << " and " << node1 << std::endl;
							std::vector<Face *> active = op.added_faces;
							flip_edges(subset, active, &bad_edges, 0);
						}
						return !bad_edges.empty();
					}


				
				// don't use edge pointer as secondary sort key, otherwise not reproducible
				struct Deterministic_sort {
						inline bool operator()(const std::pair<double, Edge *> &left, const std::pair<double, Edge *> &right) {
							return left.first < right.first;
						}
				} deterministic_sort;
				
				std::vector<Edge *> find_bad_edges(const std::vector<Edge *> &edges) {//4-1
					std::vector<std::pair<double, Edge *> > edgems;
					for (size_t e = 0; e < edges.size(); e++) {
						double m = edge_metric(edges[e]);//4-2
						if (m > 1) edgems.push_back(std::make_pair(m, edges[e]));
					}
					sort(edgems.begin(), edgems.end(), deterministic_sort);//4-1.1
					std::vector<Edge *> bad_edges(edgems.size());
					for (int e = 0; e < (int)edgems.size(); e++) bad_edges[e] = edgems[edgems.size() - e - 1].second;
					return bad_edges;
				}



				template <Space s>
				Vec3 safe_normal(Face *face) {//4-7
					if (!face) return Vec3(0);
					const Vec3 a = pos<s>(face->v[1]) - pos<s>(face->v[0]);
					const Vec3 b = pos<s>(face->v[2]) - pos<s>(face->v[0]);
					Vec3 c = cross(a, b);
					if (fabs(dot(normalize(a), normalize(b))) > 0.9 && norm(c) < 1e-6) return Vec3(0);  // unstable normal
					return normalize(c);
				}



				void project_vertex(Vert *vnew, Edge *edge, int s, double d) {//4-6
					Vert *v0 = edge_vert(edge, s, s), *v1 = edge_vert(edge, s, 1 - s);
					if (s != 0) d = 1 - d;
					Vec3 n = safe_normal<MS>(edge->adjf[0]) + safe_normal<MS>(edge->adjf[1]);//4-7
					Vec3 u = (1 - d) * v0->u + d * v1->u;
					if (norm(n) > 0.2) {
						if (!edge->n[0]->mesh->ref->raycast(u, normalize(n))) {//in 
							//shape 11/20 still undefined
							std::cout << "split:raycast failed" << std::endl;
							exit(1);
						}
					}
					vnew->u = u;
					vnew->sizing = (1.0 - d) * v0->sizing + d * v1->sizing;
				}

				void embedding_from_plasticity(const std::vector<Face *> &fs) {//4-9
					std::vector<Node *> nodes;
					std::vector<Face *> faces;
					std::vector<Edge *> edges;
					for (size_t f = 0; f < fs.size(); f++) {
						const Face *face = fs[f];
						include((Face *)face, faces);//in utol
						for (int i = 0; i < 3; i++) {
							include(face->v[i]->node, nodes);
							const Edge *edge = face->adje[i];
							include((Edge *)edge, edges);
							int s = (face == edge->adjf[0] ? 0 : 1);
							if (edge->adjf[1 - s]) {  // include adjacent "flap"
								include(edge->adjf[1 - s], faces);
								include(edge_opp_vert(edge, 1 - s)->node, nodes);
							}
						}
					}
					for (size_t n = 0; n < nodes.size(); n++) nodes[n]->y = nodes[n]->x;
					std::vector<Constraint *> no_cons;
					local_opt<PS>(nodes, faces, edges, no_cons);//4-10
				}





				void connect(Vert *vert, Node *node) {//4-8
					vert->node = node;
					include(vert, node->verts);
				}

				RemeshOp split_edge(Edge *edge, double d) { // 4-5
					Mesh &mesh = *edge->n[0]->mesh;
					RemeshOp op;
					Node *node0 = edge->n[0], *node1 = edge->n[1],
						*node = new Node((1 - d) * node0->y + d * node1->y, (1 - d) * node0->x + d * node1->x,
						(1 - d) * node0->v + d * node1->v,
							0,  // node0->label & node1->label,
							node0->flag & node1->flag, false);
					node->acceleration = (1 - d) * node0->acceleration + d * node1->acceleration;
					op.added_nodes.push_back(node);
					op.removed_edges.push_back(edge);
					op.added_edges.push_back(new Edge(node0, node, edge->theta_ideal, edge->preserve));
					op.added_edges.push_back(new Edge(node, node1, edge->theta_ideal, edge->preserve));
					Vert *vnew[2] = { NULL, NULL };
					for (int s = 0; s < 2; s++) {
						if (!edge->adjf[s]) continue;
						Vert *v0 = edge_vert(edge, s, s), *v1 = edge_vert(edge, s, 1 - s), *v2 = edge_opp_vert(edge, s);
						if (s == 0 || is_seam_or_boundary(edge)) {
							vnew[s] = new Vert(Vec3(0));
							project_vertex(vnew[s], edge, s, d);//4-6
							connect(vnew[s], node);//4-8 in mesh
							op.added_verts.push_back(vnew[s]);
						}
						else
							vnew[s] = vnew[0];
						op.added_edges.push_back(new Edge(v2->node, node, 0, 0));
						Face *f = edge->adjf[s];
						op.removed_faces.push_back(f);
						Face *nf0 = new Face(v0, vnew[s], v2, f->Sp_str, f->Sp_bend, f->material, f->damage);
						Face *nf1 = new Face(vnew[s], v1, v2, f->Sp_str, f->Sp_bend, f->material, f->damage);
						if ((std::min)(aspect(nf0), aspect(nf1)) < 1e-3) {
							op.cancel();
							return op;
						}
						op.added_faces.push_back(nf0);
						op.added_faces.push_back(nf1);
					}
					embedding_from_plasticity(op.removed_faces);//4-9
					op.apply(mesh);//5-1
					node->y = (1 - d) * node0->y + d * node1->y;
					optimize_node(node);
					plasticity_from_embedding(op.added_faces);//7-1
					local_pop_filter(op.added_faces);//8-1
					return op;
				}

				void optimize_node(Node *node) {
					std::vector<Node *> nodes(1, node);
					std::vector<Face *> faces;
					for (size_t v = 0; v < node->verts.size(); v++) append(faces, node->verts[v]->adjf);
					std::vector<Edge *> edges;
					for (size_t f = 0; f < faces.size(); f++) {
						const Face *face = faces[f];
						for (int i = 0; i < 3; i++) include(face->adje[i], edges);
					}
					std::vector<Constraint *> no_cons;
					PlasticityStash stash(faces, edges);
					local_opt<PS>(nodes, faces, edges, no_cons);
					stash.apply(faces, edges);
					local_opt<WS>(nodes, faces, edges, no_cons);
				}

				void plasticity_from_embedding(const std::vector<Face *> faces)
				{
					for (size_t n = 0; n < faces.size(); n++)
					{
						Face *face = faces[n];
						if (face->material->plastic_flow == 0)
						{
							continue;
						}
						face->Sp_str = stretch_plasticity_from_embedding(face);//7-1
					}

					std::vector<Edge*>edges;

					for (int f = 0; f < faces.size(); f++)
					{
						for (int i = 9; i < 3; i++) include(faces[f]->adje[i], edges);
					}
					for (int e = 0; e < edges.size(); e++)
					{
						Edge *edge = edges[e];
						if (!edge->adjf[0] || !edge->adjf[1]) continue;
						if (edge->adjf[0]->material->yield_curv < infinity || edge->adjf[1]->material->yield_curv < infinity)
							edge->theta_ideal = dihedral_angle<PS>(edge);
						else
							edge->theta_ideal = dihedral_angle<MS>(edge);
					}
					for (size_t f = 0; f < faces.size(); f++) recompute_Sp_bend(faces[f]);
				}


			
			/*	Mat3x3 stretch_plasticity_from_embedding(const Face *face)
				{
					double limit = face->material->plastic_limit;
					Mat3x3 F_ps = derivative(pos<PS>(face->v[0]->node), pos<PS>(face->v[1]->node), pos<PS>(face->v[2]->node), normal<PS>(face), face);
					SVD<3, 3> svd = singular_value_decomposition(F_ps);
					for (int j = 0; j < 3; j++) svd.s[j] = clamp(1.0 / svd.s[j], 1. / limit, limit);
					Mat3x3 M = svd.Vt.t() *diag(svd.s)*svd.Vt;
					return M;
				}*/


				Vert *adjacent_vert(const Node *node, const Vert *vert) {
					const Edge *edge = get_edge(node, vert->node);
					for (int i = 0; i < 2; i++)
						for (int s = 0; s < 2; s++)
							if (edge_vert(edge, s, i) == vert) return edge_vert(edge, s, 1 - i);
					return NULL;
				}


				//---------------------improve

				bool improve_some_face(MeshSubset *subset, std::vector<Face *> &active) {
					for (size_t f = 0; f < active.size(); f++) {
						Face *face = active[f];
						for (int e = 0; e < 3; e++) {
							Edge *edge = face->adje[e];
							RemeshOp op;
							if (op.empty()) op = try_edge_collapse(edge, 0);//1-10-1
							if (op.empty()) op = try_edge_collapse(edge, 1);
							if (op.empty()) continue;
							op.update(active);
							if (subset) op.update(subset->active_nodes);
							op.done();
							std::vector<Face *> fix_active = op.added_faces;
							flip_edges(subset, fix_active, 0, &active);
							return true;
						}
						remove(f--, active);
					}
					return false;
				}

				RemeshOp try_edge_collapse(Edge *edge, int which) {
					Node *node0 = edge->n[which], *node1 = edge->n[1 - which];
					if (node0->preserve || node0->label != 0 || (is_seam_or_boundary(node0) && !is_seam_or_boundary(edge)) ||
						(has_labeled_edges(node0) && !edge->preserve))
						return RemeshOp();
					if (!can_collapse(node0->mesh->parent->remeshing, edge, which)) return RemeshOp();//1-10-4
					RemeshOp op = collapse_edge(edge, which);
					if (op.empty()) return op;
					if (global_index::verbose)std::cout << "Collapsed " << node0 << " into " << node1 << std::endl;
					return op;
				}


				bool has_labeled_edges(const Node *node) {
					for (int e = 0; e < (int)node->adje.size(); e++)
						if (node->adje[e]->preserve) return true;
					return false;
				}

				//1-10-4
				bool can_collapse(Remeshing &remeshing, const Edge *edge, int i) {
					for (int s = 0; s < 2; s++) {
						const Vert *vert0 = edge_vert(edge, s, i), *vert1 = edge_vert(edge, s, 1 - i);
						if (!vert0 || (s == 1 && vert0 == edge_vert(edge, 0, i))) continue;
						for (int f = 0; f < (int)vert0->adjf.size(); f++) {
							const Face *face = vert0->adjf[f];
							if (is_in(vert1, face->v)) continue;
							const Vert *vs[3] = { face->v[0], face->v[1], face->v[2] };
							double a0 = norm(cross(vs[1]->u - vs[0]->u, vs[2]->u - vs[0]->u)) / 2;
							replace(vert0, vert1, vs);
							double a = norm(cross(vs[1]->u - vs[0]->u, vs[2]->u - vs[0]->u)) / 2;
							double asp = aspect(vs[0]->u, vs[1]->u, vs[2]->u);

							if ((a < a0 && a < 0.1 * sqr(remeshing.size_min)) || asp < remeshing.aspect_min) return false;
							for (int e = 0; e < 3; e++)
								if (vs[e] != vert1 && edge_metric(vs[NEXT(e)], vs[PREV(e)]) > 0.9) {
									return false;
								}
						}
					}
					return true;
				}



				RemeshOp collapse_edge(Edge *edge, int i) {
					/*if (is_seam_or_boundary(edge)) {
					Annotation::add(edge);
					cout << "collapse" << endl;
					cout << edge->n[i]->preserve << endl;
					wait_key();
					}*/
					Mesh &mesh = *edge->n[0]->mesh;
					RemeshOp op;
					Node *node0 = edge->n[i], *node1 = edge->n[1 - i];
					op.removed_nodes.push_back(node0);
					for (size_t e = 0; e < node0->adje.size(); e++) {
						Edge *edge1 = node0->adje[e];
						op.removed_edges.push_back(edge1);
						Node *node2 = (edge1->n[0] != node0) ? edge1->n[0] : edge1->n[1];
						if (node2 != node1 && !get_edge(node1, node2))
							op.added_edges.push_back(new Edge(node1, node2, edge1->theta_ideal, edge1->preserve));
					}
					for (int s = 0; s < 2; s++) {
						Vert *vert0 = edge_vert(edge, s, i), *vert1 = edge_vert(edge, s, 1 - i);
						if (!vert0 || (s == 1 && vert0 == edge_vert(edge, 0, i))) continue;
						op.removed_verts.push_back(vert0);
						for (size_t f = 0; f < vert0->adjf.size(); f++) {
							Face *face = vert0->adjf[f];
							op.removed_faces.push_back(face);
							if (!is_in(vert1, face->v)) {
								Vert *verts[3] = { face->v[0], face->v[1], face->v[2] };
								replace(vert0, vert1, verts);
								Face *new_face =
									new Face(verts[0], verts[1], verts[2], face->Sp_str, face->Sp_bend, face->material, face->damage);
								op.added_faces.push_back(new_face);
								// inversion test
								if (dot(normal<MS>(face), normal<MS>(new_face)) < 0) {
									op.cancel();
									return RemeshOp();
								}
								// degenerate test
								bool enforce = false;
								double asp_old = aspect(face), asp_new = aspect(new_face);
								if (asp_new < mesh.parent->remeshing.aspect_min / 4 &&
									asp_old >= mesh.parent->remeshing.aspect_min / 4 && !enforce) {
									op.cancel();
									return RemeshOp();
								}
							}
						}
					}
					// wait_key();
					embedding_from_plasticity(op.removed_faces);
					op.apply(mesh);
					plasticity_from_embedding(op.added_faces);
					local_pop_filter(op.added_faces);
					return op;
				}


				void local_pop_filter(const std::vector<Face * > &fs)
				{
					if (!global_index::magic.enable_localopt) std::cout << "local" << std::endl; return;
					std::vector<Node *> nodes;
					std::vector<Face *> faces;
					std::vector<Edge *> edges;
					for (size_t f = 0; f < fs.size(); f++)
						for (int i = 0; i < 3; i++) include(fs[f]->v[i]->node, nodes);
					for (size_t n = 0; n < nodes.size(); n++) {
						for (size_t v = 0; v < nodes[n]->verts.size(); v++) {
							const Vert *vert = nodes[n]->verts[v];
							for (size_t f = 0; f < vert->adjf.size(); f++) include(vert->adjf[f], faces);
						}
					}
					for (size_t f = 0; f < faces.size(); f++)
						for (int i = 0; i < 3; i++) include(faces[f]->adje[i], edges);
					std::vector<Constraint *> cons;
					cons = proximity_constraints(sim.cloth_meshes, sim.obstacle_meshes, sim.friction, sim.obs_friction, false, true);
					for (int h = 0; h < (int)sim.handles.size(); h++) append(cons, sim.handles[h]->get_constraints(sim.time));//in util.hpp

																															  /*if (sim.frame > 0) {
																															  for (int i=0; i<nodes.size(); i++)
																															  Annotation::add(nodes[i]);
																															  for (int i=0; i<edges.size(); i++)																												  Annotation::add(edges[i], Vec3(0,0,1));
																															  wait_key();}*/
					local_opt<WS>(nodes, faces, edges, cons);

				}


				void recompute_Sp_bend(Face *face) {
					Vec3 theta;
					for (int e = 0; e < 3; e++) theta[e] = face->adje[e]->theta_ideal;
					face->Sp_bend = edges_to_face(theta, face);  // add residual and PS reconstruction
				}





}