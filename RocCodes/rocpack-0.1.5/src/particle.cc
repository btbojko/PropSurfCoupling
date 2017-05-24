#include <cmath>

#include "particle.h"

#ifdef HAVE_OPENGL
#include "opengl.h"
#endif

std::vector<Particle*> particle;
std::vector<int> labels;

#ifdef HAVE_OPENGL
void Particle::draw() const
{
    glPushMatrix();

    glTranslatef(m_position[0], m_position[1], m_position[2]);

	float angle = deg(2.0 * acos(m_orientation[3]));

    glRotatef(angle, m_orientation[0], m_orientation[1], m_orientation[2]);

	float scale = m_time * m_growth_rate;

    glScalef(scale, scale, scale);

	/* when debugging, we highlight particles based on current event,
	 * so we do not set the color for each particle here. */
#ifndef DEBUG_DRAW
	glColor4f(r, g, b, a);
#endif

    m_shape->draw();

    glPopMatrix();
}
#endif
