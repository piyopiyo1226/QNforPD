#pragma once
#ifndef SPLINE_HPP
#define SPLINE_HPP

#include "vectors.hpp"
#include <vector>
#include <algorithm>


namespace arcsim {

	template <typename T>
	class Spline {
	public:
		// cubic Hermite spline with linear extrapolation
		struct Point {
			double t;
			T x, v;
		};
		std::vector<Point> points;
		T pos(double t) const;
		T vel(double t) const;
	};

	std::vector<double> operator+(const std::vector<double> &x, const std::vector<double> &y);
	std::vector<double> operator-(const std::vector<double> &x, const std::vector<double> &y);
	std::vector<double> operator*(const std::vector<double> &x, double a);
	std::vector<double> operator/(const std::vector<double> &x, double a);
	//
	template <typename T>
	void fill_in_velocity(Spline<T> &s, int i) {
		if (i - 1 < 0 || i + 1 >= s.points.size())
			s.points[i].v = s.points[i].x * 0.;
		else
			s.points[i].v = (s.points[i + 1].x - s.points[i - 1].x) / (s.points[i + 1].t - s.points[i - 1].t);
	}
}
#endif