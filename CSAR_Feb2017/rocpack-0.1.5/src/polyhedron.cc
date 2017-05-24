#include <cstdio>
#include <cstdlib>
#include <string.h>

#include <vector>

#include "polyhedron.h"

#ifdef HAVE_OPENGL
#include "opengl.h"
#endif

#ifdef HAVE_OPENGL
void Polyhedron::draw() const
{
	for (unsigned int i = 0; i < face.size(); i++) {
		glBegin(GL_POLYGON);
		glNormal3f(normal[i][0], normal[i][1], normal[i][2]);
		for (unsigned int j = 0; j < face[i].size(); j++) {
			const Point& v = vertex[face[i][j]];
			glVertex3f(v[0], v[1], v[2]);
		}
		glEnd();
	}
}
#endif

void Polyhedron::compute_mass_properties()
{
	/* compute face normals */

	for (unsigned int i = 0; i < face.size(); i++) {
		Vector v1 = vertex[face[i][1]] - vertex[face[i][0]];
		Vector v2 = vertex[face[i][2]] - vertex[face[i][0]];
		normal.push_back((v1 ^ v2).normalized());
		if (normal[i] * vertex[face[i][0]] < 0.0) normal[i] = -normal[i];
	}

	/* compute triangulation */
	std::vector<std::vector<int> > triangles;
	for (unsigned int n = 0; n < face.size(); n++) {
		for (unsigned int i = 1; i < face[n].size() - 1; i++) {
			std::vector<int> curr_triangle;

			curr_triangle.push_back(face[n][0]);
			curr_triangle.push_back(face[n][i]);
			curr_triangle.push_back(face[n][i+1]);

			triangles.push_back(curr_triangle);
		}
	}

	/* compute volume, center of mass and moment of inertia */
	float integral[10] = {0.0, };

	for (unsigned int n = 0; n < triangles.size(); n++) {
		float tmp0, tmp1, tmp2;
		float f1x, f2x, f3x, g0x, g1x, g2x;
		float f1y, f2y, f3y, g0y, g1y, g2y;
		float f1z, f2z, f3z, g0z, g1z, g2z;

		// Get vertices of triangle i.
		Vector& v0 = vertex[ triangles[n][0] ];
		Vector& v1 = vertex[ triangles[n][1] ];
		Vector& v2 = vertex[ triangles[n][2] ];

		// Get cross product of edges and normal vector.
		Vector V1mV0 = v1 - v0;
		Vector V2mV0 = v2 - v0;
		Vector N = V1mV0 ^ V2mV0;

		// Compute integral terms.

		tmp0 = v0[0] + v1[0];
		tmp1 = v0[0] * v0[0];
		tmp2 = tmp1 + v1[0]*tmp0;
		f1x = tmp0 + v2[0];
		f2x = tmp2 + v2[0]*f1x;
		f3x = v0[0]*tmp1 + v1[0]*tmp2 + v2[0]*f2x;
		g0x = f2x + v0[0]*(f1x + v0[0]);
		g1x = f2x + v1[0]*(f1x + v1[0]);
		g2x = f2x + v2[0]*(f1x + v2[0]);

		tmp0 = v0[1] + v1[1];
		tmp1 = v0[1]*v0[1];
		tmp2 = tmp1 + v1[1]*tmp0;
		f1y = tmp0 + v2[1];
		f2y = tmp2 + v2[1]*f1y;
		f3y = v0[1]*tmp1 + v1[1]*tmp2 + v2[1]*f2y;
		g0y = f2y + v0[1]*(f1y + v0[1]);
		g1y = f2y + v1[1]*(f1y + v1[1]);
		g2y = f2y + v2[1]*(f1y + v2[1]);

		tmp0 = v0[2] + v1[2];
		tmp1 = v0[2]*v0[2];
		tmp2 = tmp1 + v1[2]*tmp0;
		f1z = tmp0 + v2[2];
		f2z = tmp2 + v2[2]*f1z;
		f3z = v0[2]*tmp1 + v1[2]*tmp2 + v2[2]*f2z;
		g0z = f2z + v0[2]*(f1z + v0[2]);
		g1z = f2z + v1[2]*(f1z + v1[2]);
		g2z = f2z + v2[2]*(f1z + v2[2]);

		// Update integrals.
		integral[0] += N[0]*f1x;
		integral[1] += N[0]*f2x;
		integral[2] += N[1]*f2y;
		integral[3] += N[2]*f2z;
		integral[4] += N[0]*f3x;
		integral[5] += N[1]*f3y;
		integral[6] += N[2]*f3z;
		integral[7] += N[0]*(v0[1]*g0x + v1[1]*g1x + v2[1]*g2x);
		integral[8] += N[1]*(v0[2]*g0y + v1[2]*g1y + v2[2]*g2y);
		integral[9] += N[2]*(v0[0]*g0z + v1[0]*g1z + v2[0]*g2z);
	}

	integral[0] /=   6.0;
	integral[1] /=  24.0;
	integral[2] /=  24.0;
	integral[3] /=  24.0;
	integral[4] /=  60.0;
	integral[5] /=  60.0;
	integral[6] /=  60.0;
	integral[7] /= 120.0;
	integral[8] /= 120.0;
	integral[9] /= 120.0;

	/* volume  */

	m_volume = integral[0];

	/* center of mass */

	Point center = (1.0 / m_volume) * Point(integral[1], integral[2], integral[3]);

	/* change vertices to be in center of mass coordinates */

	for (unsigned int i = 0; i < vertex.size(); i++)
		vertex[i] -= center;

	/* compute bounding radius */

	radius = 0.0;
	for (unsigned int i = 0; i < vertex.size(); i++) {
		if (norm(vertex[i]) > radius)
			radius = norm(vertex[i]);
	}

#ifdef DEBUG_MASS_PROPERTIES
	printf("%s: %d: mass_properties(): radius = %f, center of mass = <%f, %f, %f> volume = %f\n",
		__FILE__, __LINE__, radius, center[0], center[1], center[2], m_volume);
#endif

#ifndef NO_SCALING
	/* scale to a bounding radius of 1.0 */

	for (unsigned int i = 0; i < vertex.size(); i++) {
		vertex[i] /= radius;
	}

	m_volume /= pow(radius, 3.0);
	radius = 1.0;
#endif

	/* moments of inertia */

	float Ixx, Ixy, Ixz, Iyx, Iyy, Iyz, Izx, Izy, Izz;

	/* in world coordinates */

	Ixx = integral[5] + integral[6];
	Ixy = -integral[7];
	Ixz = -integral[9];
	Iyx = Ixy;
	Iyy = integral[4] + integral[6];
	Iyz = -integral[8];
	Izx = Ixz;
	Izy = Iyz;
	Izz = integral[4] + integral[5];

	/* in center of mass coordinates */

	Ixx -= m_volume*(center[1]*center[1] + center[2]*center[2]);
	Ixy += m_volume*center[0]*center[1];
	Ixz += m_volume*center[2]*center[0];
	Iyx =  Ixy;
	Iyy -= m_volume*(center[2]*center[2] + center[0]*center[0]);
	Iyz += m_volume*center[1]*center[2];
	Izx =  Ixz;
	Izy =  Iyz;
	Izz -= m_volume*(center[0]*center[0] + center[1]*center[1]);

	center -= center; /* change center of mass to zero */

	m_inertia.set_value(Ixx, Ixy, Ixz, Iyx, Iyy, Iyz, Izx, Izy, Izz);

	m_inv_inertia = inverse(m_inertia);
}

void Polyhedron::write_povray_object(FILE* output) const
{
	fprintf(output, "\n#declare %s = union {\n", m_name);

	for (unsigned int i = 0; i < face.size(); i++) {
		for (unsigned int j = 1; j < face[i].size() - 1; j++) {
			const Point& a = vertex[face[i][0]];
			const Point& b = vertex[face[i][j]];
			const Point& c = vertex[face[i][j+1]];
			fprintf(output, "\ttriangle {\n");
			fprintf(output, "\t\t<% f, % f, % f>,\n",     a[0], a[1], a[2]);
			fprintf(output, "\t\t<% f, % f, % f>,\n",     b[0], b[1], b[2]);
			fprintf(output, "\t\t<% f, % f, % f>\n\t}\n", c[0], c[1], c[2]);
		}
	}
	fprintf(output, "}\n\n");
}
