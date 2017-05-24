#include "ellipsoid.h"

#ifdef HAVE_OPENGL
#include "opengl.h"
#endif

#if 0
namespace {
	static const bool registered =
		register_shape("Ellipsoid", new Ellipsoid());
}
#endif

#ifdef HAVE_OPENGL
void Ellipsoid::draw() const
{
	const int SLICES = 20; /* number of subdivisions around the z axis */
	const int STACKS = 20; /* number of subdivisions  along the z axis */

	glScalef(a, b, c);
	glutSolidSphere(1.0, SLICES, STACKS);
}
#endif
