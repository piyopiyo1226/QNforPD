#include "spline.hpp"
#include "util.hpp"


namespace arcsim
{
	template <typename T>
	static int find(const Spline<T> &s, double t) {
		int l = 0, u = s.points.size();
		while (l != u) {
			int m = (l + u) / 2;
			if (t < s.points[m].t)
				u = m;
			else
				l = m + 1;
		}
		return l;  // which is equal to u
	}

	template <typename T>
	T Spline<T>::pos(double t) const {
		int i = find(*this, t);
		if (i == 0) {
			const Point &p1 = points[i];
			return p1.x;
		}
		else if (i == points.size()) {
			const Point &p0 = points[i - 1];
			return p0.x;
		}
		else {
			const Point &p0 = points[i - 1], &p1 = points[i];
			double s = (t - p0.t) / (p1.t - p0.t), s2 = s * s, s3 = s2 * s;
			return p0.x * (2 * s3 - 3 * s2 + 1) + p1.x * (-2 * s3 + 3 * s2) +
				(p0.v * (s3 - 2 * s2 + s) + p1.v * (s3 - s2)) * (p1.t - p0.t);
		}
	}

	template <typename T>
	T Spline<T>::vel(double t) const {
		int i = find(*this, t);
		if (i == 0 || i == points.size()) {
			return T(0);
		}
		else {
			const Point &p0 = points[i - 1], &p1 = points[i];
			double s = (t - p0.t) / (p1.t - p0.t), s2 = s * s;
			return (p0.x * (6 * s2 - 6 * s) + p1.x * (-6 * s2 + 6 * s)) / (p1.t - p0.t) + p0.v * (3 * s2 - 4 * s + 1) +
				p1.v * (3 * s2 - 2 * s);
		}
	}

	std::vector<double> operator+(const std::vector<double> &x, const std::vector<double> &y) {
		std::vector<double> z((std::min)(x.size(), y.size()));
		for (int i = 0; i < (int)z.size(); i++) z[i] = x[i] + y[i];
		return z;
	}
	std::vector<double> operator-(const std::vector<double> &x, const std::vector<double> &y) {
		std::vector<double> z((std::min)(x.size(), y.size()));
		for (int i = 0; i < (int)z.size(); i++) z[i] = x[i] - y[i];
		return z;
	}
	std::vector<double> operator*(const std::vector<double> &x, double a) {
		std::vector<double> y(x.size());
		for (int i = 0; i < (int)y.size(); i++) y[i] = x[i] * a;
		return y;
	}
	std::vector<double> operator/(const std::vector<double> &x, double a) { return x * (1 / a); }

	template class Spline<Vec3>;
	template class Spline<Transformation>;
	template class Spline<double>;
	template class Spline<std::vector<double> >;



}