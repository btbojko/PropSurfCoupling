#include "sphere.h"

#ifdef HAVE_OPENGL
#include "opengl.h"
#endif

#if 0
namespace {
	static const bool registered =
		register_shape("sphere", new Sphere());
}
#endif

#ifdef HAVE_OPENGL
void Sphere::draw() const
{
	const int SLICES = 20; /* number of subdivisions around the z axis */
	const int STACKS = 20; /* number of subdivisions  along the z axis */

	glutSolidSphere(radius, SLICES, STACKS);
}
#endif
