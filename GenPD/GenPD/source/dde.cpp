#include "dde.hpp"
#include <algorithm>
#include <math.h>

namespace arcsim {
	Vec4 stretching_stiffness(const Mat2x2 &G, const StretchingSamples &samples) {
		double a = (G(0, 0) + 0.25) * nsamples;
		double b = (G(1, 1) + 0.25) * nsamples;
		double c = std::abs(G(0, 1)) * nsamples;
		a = clamp(a, 0.0, nsamples - 1 - 1e-5);
		b = clamp(b, 0.0, nsamples - 1 - 1e-5);
		c = clamp(c, 0.0, nsamples - 1 - 1e-5);
		int ai = (int)floor(a);
		int bi = (int)floor(b);
		int ci = (int)floor(c);
		if (ai < 0) ai = 0;
		if (bi < 0) bi = 0;
		if (ci < 0) ci = 0;
		if (ai > nsamples - 2) ai = nsamples - 2;
		if (bi > nsamples - 2) bi = nsamples - 2;
		if (ci > nsamples - 2) ci = nsamples - 2;
		a = a - ai;
		b = b - bi;
		c = c - ci;
		double weight[2][2][2];
		weight[0][0][0] = (1 - a) * (1 - b) * (1 - c);
		weight[0][0][1] = (1 - a) * (1 - b) * (c);
		weight[0][1][0] = (1 - a) * (b) * (1 - c);
		weight[0][1][1] = (1 - a) * (b) * (c);
		weight[1][0][0] = (a) * (1 - b) * (1 - c);
		weight[1][0][1] = (a) * (1 - b) * (c);
		weight[1][1][0] = (a) * (b) * (1 - c);
		weight[1][1][1] = (a) * (b) * (c);
		Vec4 stiffness = Vec4(0);
		for (int i = 0; i < 2; i++)
			for (int j = 0; j < 2; j++)
				for (int k = 0; k < 2; k++)
					for (int l = 0; l < 4; l++) {
						stiffness[l] += samples.s[ai + i][bi + j][ci + k][l] * weight[i][j][k];
					}
		return stiffness;
	}

	double bending_stiffness(const Edge *edge, int side, const BendingData &data, double l, double theta,
		double initial_angle) {
		double curv = theta * l / (edge->adjf[0]->a + edge->adjf[1]->a);
		double alpha = curv / 2;
		double value = alpha * 0.2;  // because samples are per 0.05 cm^-1 = 5 m^-1
		if (value > 4) value = 4;
		int value_i = (int)value;
		if (value_i < 0) value_i = 0;
		if (value_i > 3) value_i = 3;
		value -= value_i;
		Vec3 du = edge_vert(edge, side, 1)->u - edge_vert(edge, side, 0)->u;
		double bias_angle = (atan2f(du[1], du[0]) + initial_angle) * 4 / M_PI;
		if (bias_angle < 0) bias_angle = -bias_angle;
		if (bias_angle > 4) bias_angle = 8 - bias_angle;
		if (bias_angle > 2) bias_angle = 4 - bias_angle;
		int bias_id = (int)bias_angle;
		if (bias_id < 0) bias_id = 0;
		if (bias_id > 1) bias_id = 1;
		bias_angle -= bias_id;
		double actual_ke = data.d[bias_id][value_i] * (1 - bias_angle) * (1 - value) +
			data.d[bias_id + 1][value_i] * (bias_angle) * (1 - value) +
			data.d[bias_id][value_i + 1] * (1 - bias_angle) * (value)+
			data.d[bias_id + 1][value_i + 1] * (bias_angle) * (value);
		if (actual_ke < 0) actual_ke = 0;
		return actual_ke;
	}
}
