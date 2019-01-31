#pragma once

#include "geometry.hpp"

namespace arcsim {

	static const int nsamples = 30;



	typedef Vec<4> Vec4;
	struct StretchingSamples {
		Vec4 s[40][40][40];
	};

	struct BendingData {
		double d[3][5];
	};

	Vec4 stretching_stiffness(const Mat2x2 &G, const StretchingSamples &samples);
	double bending_stiffness(const Edge *edge, int side, const BendingData &data, double l, double theta,
		double initial_angle = 0);


}