#ifndef _SPHERE_H
#define _SPHERE_H

#include "point.h"
#include "matrix.h"
#include "convex.h"

class Sphere : public Convex {
  public:
	Sphere() : radius(1.0) {}
	Sphere(float radius) : radius(radius) {}
	~Sphere() {}

	ShapeType type() const { return SPHERE; }
	const char* name() const { return "sphere"; }

	float bounding_radius() const { return radius; }

	Point support(const Vector& v) const {
		float s = norm(v);
		if (s > EPSILON)
			return v * (radius / s);
		else
			return Point(0.0, 0.0, 0.0);
	}

	Matrix inertia() const {
		return Matrix(0.4 * radius*radius);
	}

	Matrix inverse_inertia() const {
		return Matrix(2.5 / (radius*radius));
	}

	float volume() const {
		return 4.0/3.0 * M_PI * radius*radius*radius;
	}

#ifdef HAVE_OPENGL
	void draw() const;
#endif

  private:
	float radius;
};

#endif
