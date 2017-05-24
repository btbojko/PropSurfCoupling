#include "cylinder.h"

#ifdef HAVE_OPENGL
#include "opengl.h"
#endif

#ifdef HAVE_OPENGL
void Cylinder::draw() const
{
	static GLUquadricObj* cylinder = gluNewQuadric();

	gluQuadricDrawStyle(cylinder, GLU_FILL);

	glPushMatrix();
	glBegin(GL_TRIANGLE_FAN);
	glVertex3f(0.0, h/2.0, 0.0);
	for (int i = 0; i <= 50; i++) {
		float t = i/50.0 * 2 * M_PI;
		glVertex3f(r * sin(t), h/2.0, r * cos(t));
	}
	glEnd();
	glBegin(GL_TRIANGLE_FAN);
	glNormal3f(0.0, -1.0, 0.0);
	glVertex3f(0.0, -h/2.0, 0.0);
	for (int i = 0; i <= 50; i++) {
		float t = i/50.0 * 2 * M_PI;
		glVertex3f(r * sin(t), -h/2.0, r * cos(t));
	}
	glEnd();
	glRotated(90.0, 1.0, 0.0, 0.0);
	glTranslated(0.0, 0.0, -h/2.0);
	gluCylinder(cylinder, r, r, h, 50, 50);
	glPopMatrix();
}
#endif
