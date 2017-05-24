#ifndef _ELLIPSOID_H
#define _ELLIPSOID_H

#include <cstring>

#include "point.h"
#include "matrix.h"
#include "convex.h"

class Ellipsoid : public Convex {
  public:
	Ellipsoid() : a(1.0), b(1.0), c(1.0), r(1.0)
	{
		const char *default_name = "e_sphere";

		m_name = default_name;
	}

	Ellipsoid(const char *name, float a, float b, float c) : a(a), b(b), c(c)
	{
		m_name = strdup(name); r = max<float>(a, b, c);
		a /= r; b /= r; c /= r; r = 1.0;
	}

	~Ellipsoid() {}

	ShapeType type() const { return CONVEX; }
	const char* name() const { return m_name; }

	float bounding_radius() const { return r; }

	Point support(const Vector& v) const {
		if (norm(v) > EPSILON) {
			Vector q(a*a*v[0], b*b*v[1], c*c*v[2]);
			return (1.0 / sqrt(v*q)) * q;
		} else
			return Point(0.0, 0.0, 0.0);
	}

	Matrix inertia() const {
		return 0.2 * Matrix(b*b + c*c, a*a + c*c, a*a + b*b);
	}

	Matrix inverse_inertia() const {
		return 5.0 * Matrix(1.0/(b*b + c*c), 1.0/(a*a + c*c), 1.0/(a*a + b*b));
	}

	float volume() const {
		return 4.0/3.0 * M_PI * a * b * c;
	}

#ifdef HAVE_OPENGL
	void draw() const;
#endif

  private:
	float a, b, c, r; const char *m_name;
};

#endif
