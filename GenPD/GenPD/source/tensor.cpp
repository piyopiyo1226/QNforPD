#include "tensor.hpp"

namespace arcsim
{
	template <typename T>
	T head(const std::vector<T> &v) {
		return v.front();
	}
	template <typename T>
	std::vector<T> tail(const std::vector<T> &v) {
		std::vector<T> w(v.size() - 1);
		for (int i = 0; i < w.size(); i++) w[i] = v[i + 1];
		return w;
	}
	template <typename T>
	std::vector<T> cons(const T &x, const std::vector<T> &v) {
		std::vector<T> w(v.size() + 1);
		w[0] = x;
		for (int i = 1; i < w.size(); i++) w[i] = v[i - 1];
		return w;
	}
	struct Disk {
		Vec2 c;
		double r;
		Disk()
			: c(Vec2(0))
			, r(0) {}
		Disk(const Vec2 &c, double r)
			: c(c)
			, r(r) {}
	};

	bool enclosed(const Disk &disk0, const Disk &disk1) { return norm(disk0.c - disk1.c) + disk0.r <= disk1.r + 1e-6; }

	Disk minidisk(const std::vector<Disk> &P);
	Disk b_minidisk(const std::vector<Disk> &P, const std::vector<Disk> &R);
	Disk b_md(const std::vector<Disk> &R);
	Disk apollonius(const Disk &disk1, const Disk &disk2, const Disk &disk3);
	Disk welzls_algorithm(const std::vector<Disk> &disks) { return minidisk(disks); }


	Disk b_minidisk(const std::vector<Disk> &P, const std::vector<Disk> &R) {
		if (P.empty() || R.size() == 3) return b_md(R);
		Disk p = head(P);
		std::vector<Disk> P_ = tail(P);
		Disk D = b_minidisk(P_, R);
		if (enclosed(p, D))
			return D;
		else
			return b_minidisk(P_, cons(p, R));
	}
	Disk minidisk(const std::vector<Disk> &P) {
		if (P.empty()) return Disk();
		Disk p = head(P);
		std::vector<Disk> P_ = tail(P);
		Disk D = minidisk(P_);
		if (enclosed(p, D))
			return D;
		else
			return b_minidisk(P_, std::vector<Disk>(1, p));
	}
	Disk b_md(const std::vector<Disk> &R) {
		if (R.empty())
			return Disk();
		else if (R.size() == 1)
			return head(R);
		else if (R.size() == 2) {
			double d = norm(R[0].c - R[1].c);
			double r = (R[0].r + d + R[1].r) / 2;
			double t = (r - R[0].r) / d;
			return Disk(R[0].c + t * (R[1].c - R[0].c), r);
		}
		else
			return apollonius(R[0], R[1], R[2]);
	}

	Disk apollonius(const Disk &disk1, const Disk &disk2, const Disk &disk3) {
		// nicked from http://rosettacode.org/mw/index.php?title=Problem_of_Apollonius&oldid=88212
#define DEFXYR(N)               \
    double x##N = disk##N.c[0]; \
    double y##N = disk##N.c[1]; \
    double r##N = disk##N.r;
		DEFXYR(1);
		DEFXYR(2);
		DEFXYR(3);
#undef DEFXYR
		int s1 = 1, s2 = 1, s3 = 1;
		double v11 = 2 * x2 - 2 * x1;
		double v12 = 2 * y2 - 2 * y1;
		double v13 = x1 * x1 - x2 * x2 + y1 * y1 - y2 * y2 - r1 * r1 + r2 * r2;
		double v14 = 2 * s2 * r2 - 2 * s1 * r1;
		double v21 = 2 * x3 - 2 * x2;
		double v22 = 2 * y3 - 2 * y2;
		double v23 = x2 * x2 - x3 * x3 + y2 * y2 - y3 * y3 - r2 * r2 + r3 * r3;
		double v24 = 2 * s3 * r3 - 2 * s2 * r2;
		double w12 = v12 / v11;
		double w13 = v13 / v11;
		double w14 = v14 / v11;
		double w22 = v22 / v21 - w12;
		double w23 = v23 / v21 - w13;
		double w24 = v24 / v21 - w14;
		double P = -w23 / w22;
		double Q = w24 / w22;
		double M = -w12 * P - w13;
		double N = w14 - w12 * Q;
		double a = N * N + Q * Q - 1;
		double b = 2 * M * N - 2 * N * x1 + 2 * P * Q - 2 * Q * y1 + 2 * s1 * r1;
		double c = x1 * x1 + M * M - 2 * M * x1 + P * P + y1 * y1 - 2 * P * y1 - r1 * r1;
		double D = b * b - 4 * a * c;
		double rs = (-b - sqrt(D)) / (2 * a);
		double xs = M + N * rs;
		double ys = P + Q * rs;
		return Disk(Vec2(xs, ys), rs);
	}

	Mat2x2 tensor_max(const std::vector<Mat2x2> &Ms) {
		int n = Ms.size();
		std::vector<Disk> disks;
		for (int i = 0; i < n; i++) {
			const Mat2x2 &M = Ms[i];
			if (trace(M) == 0) continue;
			disks.push_back(Disk(Vec2((M(0, 0) - M(1, 1)) / 2, (M(0, 1) + M(1, 0)) / 2), (M(0, 0) + M(1, 1)) / 2));
		}
		Disk disk = welzls_algorithm(disks);
		return disk.c[0] * Mat2x2(Vec2(1, 0), Vec2(0, -1)) + disk.c[1] * Mat2x2(Vec2(0, 1), Vec2(1, 0)) +
			disk.r * Mat2x2(Vec2(1, 0), Vec2(0, 1));
	}

}