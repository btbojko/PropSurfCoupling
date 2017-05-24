#ifndef _CYLINDER_H
#define _CYLINDER_H

#include <cstring>

#include "point.h"
#include "matrix.h"
#include "convex.h"

class Cylinder : public Convex {
  public:
	Cylinder() : r(1.0), h(1.0)
	{
		const char *default_name = "cylinder";

		m_name = default_name;
	}

	Cylinder(const char *name, float r, float h) : r(r), h(h)
	{
		m_name = strdup(name);
	}

	~Cylinder() {}

	ShapeType type() const { return CONVEX; }
	const char* name() const { return m_name; }

	float bounding_radius() const { return sqrt(r*r+h*h/4.0); }

	Point support(const Vector& v) const {
		float rho = sqrt(v[0]*v[0] + v[2]*v[2]);
		float   y = v[1] < 0 ? -h/2.0 : h/2.0;
		return Vector(r*v[0]/rho, y, r*v[2]/rho);
	}

	Matrix inertia() const {
		return Matrix(1.0/12.0 * (3*r*r+h*h), 1.0/12.0 * (3*r*r+h*h), 0.5*r*r);
	}

	Matrix inverse_inertia() const {
		return Matrix(12.0 / (3*r*r+h*h),12.0 / (3*r*r+h*h),2.0/(r*r));
	}

	float volume() const {
		return M_PI * r * r * h;
	}

#ifdef HAVE_OPENGL
	void draw() const;
#endif

  private:
	float r, h; const char *m_name;
};

#endif
