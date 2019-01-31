#pragma once
#include "vectors.hpp"
#include "spline.hpp"


namespace arcsim
{
	// Transform the mesh
	struct Quaternion {
		double s;
		Vec3 v;
		Vec3 rotate(const Vec3 &point) const;
		static Quaternion from_axisangle(const Vec3 &axis, double angle);
		std::pair<Vec3, double> to_axisangle() const;
		Quaternion operator+(const Quaternion &q) const;
		Quaternion operator-(const Quaternion &q) const;
		Quaternion operator-() const;
		Quaternion operator*(const Quaternion &q) const;
		Quaternion operator*(double scalar) const;
		Quaternion operator/(double scalar) const;
	};

	Quaternion normalize(const Quaternion &q);
	Quaternion inverse(const Quaternion &q);
	double norm2(const Quaternion &q);
	inline std::ostream &operator<<(std::ostream &out, const Quaternion &q) {
		out << "(" << q.s << ", " << q.v << ")";
		return out;
	}

	struct Transformation {
		Vec3 translation;
		double scale;
		Quaternion rotation;
		Transformation(double factor = 1);
		Vec3 apply(const Vec3 &point) const;
		Vec3 apply_vec(const Vec3 &vec) const;
		Transformation operator+(const Transformation &t) const;
		Transformation operator-(const Transformation &t) const;
		Transformation operator*(const Transformation &t) const;
		Transformation operator*(double scalar) const;
		Transformation operator/(double scalar) const;
	};

	Transformation identity();
	//Transformation identity() { return Transformation(); }

	Transformation inverse(const Transformation &tr);
	inline std::ostream &operator<<(std::ostream &out, const Transformation &t) {
		out << "(translation: " << t.translation << ", rotation: " << t.rotation << ", scale: " << t.scale << ")";
		return out;
	}

	typedef Spline<Transformation> Motion;
	typedef std::pair<Transformation, Transformation> DTransformation;

	void clean_up_quaternions(Motion &motion);  // remove sign flips

	Transformation get_trans(const Motion &motion, double t);
	DTransformation get_dtrans(const Motion &motion, double t);
	Vec3 apply_dtrans(const DTransformation &dT, const Vec3 &x0, Vec3 *vel = NULL);
	Vec3 apply_dtrans_vec(const DTransformation &dT, const Vec3 &v0);

}