

////---------------------end of structs-----------------------------------
//
//
//
//
//
//struct Plane//for collision
//{
//	Plane() {}
//	Plane(const Vec3 &x0, const Vec3 &n):x0(x0),n(n){}
//	Vec3 x0, n;
//};
//
//class CollisionProxy
//{
//public: 
//	virtual ~CollisionProxy() {  };
//	virtual CollisionProxy *clone(Mesh &mesh) = 0; //Ç±Ç±ÇÃèëÇ´ï˚íçà”
//	virtual void update(Mesh &mesh) = 0;
//	virtual Constraint *constraint(const Node *node) = 0;
//};
//
//
//const double infinity = std::numeric_limits<double>::infinity();
//
//class FloorProxy : public CollisionProxy
//{
//public:
//	FloorProxy(Mesh &mesh) {	update(mesh);}
//	Constraint *constraint(const Node *node)
//	{
//		if (!is_free(node) || node->x[1] - center.x[1] > ::magic.repulsion_thickness) return 0;
//	
//		IneqCon *con = new IneqCon;
//		con->nodes[0] = (Node *)node;
//		con->nodes[1] = &center;
//		con->nodes[2] = 0;
//		con->nodes[3] = 0;
//
//		for (int n = 0; n < 4; n++) con->free[n] = n == 0;
//
//		con->w[0] = 1;
//		con->w[1] = -1;
//		con->w[2] = 0;
//		con->w[3] = 0;
//
//		con->stiff = ::magic.collision_stiffness* node->a;
//		con->n = Vec3(0, 1, 0);
//		//con->mu = sim.obs_friction; temp
//
//		return con;
//
//	};
//
//	CollisionProxy *clone(Mesh &mesh)
//	{
//		return new FloorProxy(mesh);
//	}
//	void update(Mesh &mesh)
//	{
//		double y = infinity;
//		for (size_t n = 0; n < mesh.nodes.size(); n++)
//			if (mesh.nodes[n]->x[1] < y)y = mesh.nodes[n]->x[1];//calc for min
//		center.x = Vec3(0.0, y, 0.0);
//
//	};
//private:
//	Node center;
//};
//
//
//
//
////for add the collision 
//
//
//
//

//
//	
//
//
//	
//
//
//
////dynamic remesh------------------------------------------------------------------------------------//
////------------------------------------------------------------------
//
//
//

//
//	void RemeshOp::update(std::vector<Face *> &v) {
//		exclude_all(removed_faces, v);
//		include_all(added_faces, v);
//	}
//
//	void RemeshOp::update(std::vector<Edge *> &v) {
//		exclude_all(removed_edges, v);
//		include_all(added_edges, v);
//	}
//
//	void RemeshOp::update(std::vector<Node *> &v) {
//		exclude_all(removed_nodes, v);
//		include_all(added_nodes, v);
//	}
//

////------------------------------------for dynamic remesh-------------------------------------------//
//
////
//
//

//

//
//	//-----------------------------------------3--main flip -------------------------------------------------//
//	
//


//

//
//
//	//void flip_edges(MeshSubset *subset, std::vector<Face *> &active_faces, std::vector<Edge *> *update_edges,
//	//	std::vector<Face *> *update_faces) {
//	//	int N = 3 * active_faces.size();
//	//	for (int i = 0; i < N; i++) {  // don't loop without bound
//	//		if (!flip_some_edges(subset, active_faces, update_edges, update_faces)) return;
//	//	}
//	//}
//
//	//----------------------------------------------3 end----------------------------------------------//
//
//	//-------------------------------------------------4start ----------------------------------------------//
//
//
//

//
//
////4-1------------------------------
////	split_worst_edge main for 4
//

////
//	//-------------------------------------------------4 end--------------------------//
//
//
////--------------------------------11------------------------------------------//



//
////-------------------------------------------------------//
//
//
//
////
//////---------------------------------------separate---------------------
////struct Ixn {  // intersection
////	Face *f0, *f1;
////	double l;
////	Vec3 g0[3], g1[3];  // dl/dx for each of the faces' nodes
////	Ixn() {}
////	Ixn(const Face *f0, const Face *f1)
////		: f0((Face *)f0)
////		, f1((Face *)f1) {}
////};
////
//////8-1-3
////
////
////std::vector<Ixn> find_intersections(const std::vector<AccelStruct *> &accs, const std::vector<AccelStruct *> &obs_accs) {
////	/*if (!::ixns) {
////		::nthreads = omp_get_max_threads();
////		::ixns = new std::vector<Ixn>[::nthreads];
////	}
////	for (int t = 0; t < ::nthreads; t++) ::ixns[t].clear();
////	for_overlapping_faces(accs, obs_accs, ::thickness, find_face_intersection);
////	std::vector<Ixn> ixns;
////	for (int t = 0; t < ::nthreads; t++) append(ixns, ::ixns[t]);
////	return ixns;*/
////}
////
////void update_x0(Mesh &mesh) {
////	for (int n = 0; n < (int)mesh.nodes.size(); n++) mesh.nodes[n]->x0 = mesh.nodes[n]->x;
////}
////void separate(std::vector<Mesh *> &meshes, const std::vector<Mesh *> &obs_meshes) {
////	if (false) {
////		// Unit tests.
////		Vec3 e0, e1;
////		Vec3 f0, f1, f2, fn;
////		Vec3 pt;
////
////		f0 = Vec3(1, 1, 0);
////		f1 = Vec3(-1, 1, 0);
////		f2 = Vec3(0, -1, 0);
////		fn = Vec3(0, 0, 1);
////
////		e0 = Vec3(0, 0.5, -0.5);
////		e1 = Vec3(0, 0.5, 0.5);
////	//	assert(edge_face_intersection(e0, e1, f0, f1, f2, fn, pt));
////		e0 = Vec3(1, 1, -1.0);
////		e1 = Vec3(1, 1, 1.0);
////	//	assert(edge_face_intersection(e0, e1, f0, f1, f2, fn, pt));
////		e0 = Vec3(0, 0.5, 0.1);
////		e1 = Vec3(0, 0.5, 1.0);
////	//	assert(!edge_face_intersection(e0, e1, f0, f1, f2, fn, pt));
////		e0 = Vec3(1, 0.5, -0.5);
////		e1 = Vec3(1, 0.5, 0.5);
////		//assert(!edge_face_intersection(e0, e1, f0, f1, f2, fn, pt));
////		// End of unit tests.
////	}
////
////	std::vector<AccelStruct *> accs = create_accel_structs(meshes, false),
////		obs_accs = create_accel_structs(obs_meshes, false);
////
////	// double D = ::magic.separation_step_size;  // step size
////	// //std::cout << "Separation with step_size: " << D << std::endl;
////
////	std::map<Node *, Vec3> g;
////
////	::obs_mass = 1e3;
////	for (int deform = 0; deform <= 1; deform++) {
////		::deform_obstacles = deform;
////
////		for (size_t iter = 0; iter < 100; ++iter) {
////			const bool print_logs = true;  // (iter % 16 == 0);
////
////			std::vector<Ixn> ixns = find_intersections(accs, obs_accs);
////			if (ixns.empty()) {
////				break;
////			}
////
////			// Aggregate the gradient vectors.
////			double l_total = 0;
////			g.clear();
////			for (size_t i = 0; i < ixns.size(); i++) {
////				Ixn &ixn = ixns[i];
////				l_total += ixn.l;
////				for (int v = 0; v < 3; v++) {
////					g[ixn.f0->v[v]->node] += ixn.g0[v];
////					g[ixn.f1->v[v]->node] += ixn.g1[v];
////				}
////			}
////
////			if (l_total == 0) break;
////
////			// for (size_t i = 0; i < ixns.size(); i++) {
////			//     if (ixns[i].l == 0)
////			//         continue;
////			//     Annotation::add(ixns[i].f0);
////			//     Annotation::add(ixns[i].f1);
////			// }
////			// wait_key();
////
////			if (print_logs)
////				std::cout << "iter: " << iter << " n_ixns: " << ixns.size() << " l_total: " << l_total << std::endl;
////
////			std::vector<std::vector<Node *> > components = connected_components(ixns);
////
////			// for (size_t c = 0; c < components.size(); c++) {
////			//     double t = 2*M_PI*c/components.size();
////			//     Vec3 color = (Vec3(cos(t),cos(t+2*M_PI/3),cos(t+4*M_PI/3))
////			//                   + Vec3(1))/2.;
////			//     for (size_t i = 0; i < components[c].size(); i++)
////			//         Annotation::add(components[c][i], color);
////			// }
////			// wait_key();
////
////			// Find per-component displacement x = C y
////			// which minimizes 1/2 x' M x = 1/2 y' Mc y
////			// such that -l >= g' x = gc' y
////			// So y = -l (Mc^-1 gc) / (gc' Mc^-1 gc)
////
////			// Compute inverse-mass-weighted gradient
////			std::vector<double> mc(components.size());  // mass of component
////			std::vector<Vec3> gc(components.size());    // per-component gradient
////			double denom = 0;
////			for (size_t c = 0; c < components.size(); c++) {
////				const std::vector<Node *> &component = components[c];
////				mc[c] = 0;
////				gc[c] = Vec3(0);
////				for (size_t n = 0; n < component.size(); n++) {
////					Node *node = component[n];
////					if (!is_free(node) && !deform_obstacles) continue;
////					mc[c] += get_mass(node);
////					gc[c] += g[node];
////				}
////				denom += norm2(gc[c]) / mc[c];
////			}
////			double step_length = (l_total + ::thickness) / denom;
////
////			// Apply the displacements.
////			for (size_t c = 0; c < components.size(); c++) {
////				const std::vector<Node *> &component = components[c];
////				Vec3 displacement = -step_length * gc[c] / mc[c];
////				for (size_t n = 0; n < component.size(); n++) {
////					Node *node = component[n];
////					if (!is_free(node) && !deform_obstacles) continue;
////					node->x += displacement;
////				}
////			}
////
////			// Update the world-space information.
////			for (size_t m = 0; m < meshes.size(); m++) {
////				compute_ws_data(*meshes[m]);
////				update_accel_struct(*accs[m]);
////			}
////			for (size_t o = 0; o < obs_meshes.size(); o++) compute_ws_data(*obs_meshes[o]);
////			for (size_t o = 0; o < obs_accs.size(); o++) update_accel_struct(*obs_accs[o]);
////		}
////
////		if (deform_obstacles) ::obs_mass /= 2;
////	}
////
////	// Wrap-up.
////	for (int m = 0; m < (int)meshes.size(); m++) {
////		// Save collision-free configurations.
////		update_x0(*meshes[m]);
////	}
////	for (int o = 0; o < (int)obs_meshes.size(); o++) {
////		// Save collision-free configurations.
////		update_x0(*obs_meshes[o]);
////	}
////	destroy_accel_structs(accs);
////	destroy_accel_structs(obs_accs);
////	// std::cout << "Separation done." << std::endl;
////}
////
////
////
////
//////----------------------------------------separate end---------------------------------------------------------------//
//////--------------------global ä÷êî
////
////
//
//
