#ifndef _POLYHEDRON_H
#define _POLYHEDRON_H

#include <vector>
#include <algorithm>
#include <cstdio>
#include <cstring>

#include "point.h"
#include "matrix.h"
#include "convex.h"

class Polyhedron : public Convex {
  public:
	Polyhedron(const char *name,
	           std::vector<Point> vertices,
			   std::vector<std::vector<int> > faces)
		: vertex(vertices), face(faces)
	{
		m_name = strdup(name); compute_mass_properties();
	}

	~Polyhedron() { free(m_name); vertex.clear(); face.clear(); }

	ShapeType type() const { return CONVEX; }
	const char* name() const { return m_name; }
	float bounding_radius() const { return radius; }

	Point support(const Vector& v) const
	{
		int c = 0;
		float h = v * vertex[0];

		for (unsigned int i = 1; i < vertex.size(); i++) {
			float d = v * vertex[i];
			if (d > h) {
				h = d; c = i;
			}
		}

		return vertex[c];
	}

	Matrix inertia() const
	{
		return m_inertia;
	}

	Matrix inverse_inertia() const
	{
		return m_inv_inertia;
	}

	float volume() const
	{
		return m_volume;
	}

	void write_povray_object(FILE* output) const;

#ifdef HAVE_OPENGL
	void draw() const;
#endif

  private:
	float radius;
	float m_volume;
	std::vector<Point> vertex;
	std::vector<Vector> normal;
	std::vector<std::vector<int> > face;
	Matrix m_inertia, m_inv_inertia;
	char *m_name;

	void compute_mass_properties();
};

#endif
