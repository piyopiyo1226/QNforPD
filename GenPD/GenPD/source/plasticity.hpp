#pragma once
#include "mesh_struct.hpp"

namespace arcsim {
	Mat3x3 edges_to_face(const Vec3 &theta, const Face *face);
	Mat3x3 stretch_plasticity_from_embedding(const Face *face);

}