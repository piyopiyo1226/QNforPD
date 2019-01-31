//#include "obstacle.hpp"
//
//Mesh &Obstacle::get_mesh() { return curr_state_mesh; }
//
//const Mesh &Obstacle::get_mesh() const { return curr_state_mesh; }
//
//Mesh &Obstacle::get_mesh(double time) {
//	if (time > end_time) delete_mesh(curr_state_mesh);
//	if (time < start_time || time > end_time) return curr_state_mesh;
//	if (!activated) curr_state_mesh = deep_copy(base_mesh);
//	if (transform_spline) {
//		DTransformation dtrans = get_dtrans(*transform_spline, time);
//		Mesh &mesh = curr_state_mesh;
//		for (int n = 0; n < (int)curr_state_mesh.nodes.size(); n++)
//			mesh.nodes[n]->x = apply_dtrans(dtrans, base_mesh.nodes[n]->x, &mesh.nodes[n]->v);
//		compute_ws_data(mesh);
//		if (mesh.proxy) mesh.proxy->update(mesh);
//	}
//	if (!activated) update_x0(curr_state_mesh);
//	activated = true;
//	for (int i = 0; i < curr_state_mesh.nodes.size(); i++) curr_state_mesh.nodes[i]->mesh = &curr_state_mesh;
//	return curr_state_mesh;
//}
//
