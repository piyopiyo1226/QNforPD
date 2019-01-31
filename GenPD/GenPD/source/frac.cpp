#pragma once

#include "frac.hpp"
#include "remesh.hpp"
#include "physics.hpp"

namespace arcsim{
Mat3x3 compute_sigma(const Face *face) { /*
										 Mat3x3 F_str = deformation_gradient<WS>(face);
										 Mat3x3 G_str = (F_str.t()*F_str - Mat3x3(1)) * 0.5;
										 Mat3x3 sigma_str = material_model(face, G_str);

										 // add explicit bending
										 Mat3x3 F_bend = Mat3x3(1);
										 for (int i=0; i<3; i++)
										 F_bend += face->v[i]->node->curvature * (0.5 * face->material->fracture_bend_thickness);
										 Mat3x3 G_bend = (F_bend.t()*F_bend - Mat3x3(1)) * 0.5;
										 Mat3x3 sigma_bend = material_model(face, G_bend);

										 return get_positive(sigma_str) + get_positive(sigma_bend);*/
	Mat3x3 F_str = deformation_gradient<WS>(face);
	// Mat3x3 G_str = (F_str.t()*F_str - Mat3x3(1)) * 0.5;
	// Mat3x3 sigma_str = material_model(face, G_str);

	// add explicit bending
	Mat3x3 F_bend = Mat3x3(1);
	for (int i = 0; i < 3; i++) F_bend += face->v[i]->node->curvature * (0.5 * face->material->fracture_bend_thickness);
	F_bend = F_str * F_bend;
	Mat3x3 G_bend = (F_bend.t() * F_bend - Mat3x3(1)) * 0.5;
	Mat3x3 sigma_bend = material_model(face, G_bend);

	return get_positive(sigma_bend);
}
}