#include <cstdio>

#include <vector>
#include <utility>

#include "stime.h"
#include "event.h"
#include "hgrid.h"
#include "particle.h"
#include "boundary.h"
#include "collision.h"

#include "gjk.h"

Event::Event(unsigned int id, Particle* A)
	: m_id(id), A(A), cA(A->collisions())
{
	predict();
}

unsigned int Event::secondary_id() const
{
	return B ? B->get_event_id() : -1;
}

void Event::print()
{
	static char types[] = { 'I', 'C', 'B', 'T' };
	fprintf(stderr, "%c %.20f ", types[m_type], m_time);
	fprintf(stderr, "%d%c", id(), cA == A->collisions() ? ' ' : '*');
	if (B)
		fprintf(stderr, "%d%c\n", secondary_id(), cB == B->collisions() ? ' ' : '*');
	else
		fprintf(stderr, "\n");
}

void Event::process()
{
	bool validA = A->collisions() == cA;
	bool validB = B && B->collisions() == cB;

	switch (m_type) {

		case COLLISION: {
			A->sync(m_time);
			B->sync(m_time);

			if (validA && validB) {
				if (boundary->is_periodic())
					A->translate(m_shift);

				A->collided(); B->collided();

				collision_response(*A, *B);

				if (boundary->is_periodic()) {
					((PeriodicBoundary*)(boundary))->reposition(*A);
				}

			}
			break;
		}

		case BOUNDARY_COLLISION: {
			A->sync(m_time);
			if (validA) {
				boundary->collision_response(*A);
				A->collided();
			}
			break;
		}

		case TRANSFER: {
			A->sync(m_time);
			hgrid->rehash(A);
			break;
		}

		default:
			return;
	}
}

void Event::predict()
{
	static std::vector<Particle*> neighbors;

	/* particle A and id never change */
	cA = A->collisions(); B = NULL; cB = 0; m_shift = Vector(0.0, 0.0, 0.0);

	/* check collision with boundary */
	m_type = BOUNDARY_COLLISION;
	m_time = boundary->predict_collision(*A);

	/* check cell transfer time */
	float t_transfer = hgrid->next_rehash_time(A);

	if (t_transfer < m_time) {
		m_type = TRANSFER;
		m_time = t_transfer;
	}

	/* check collision for particle and its reflections if periodic */
	if (!boundary->is_periodic()) {
		neighbors.clear();
		hgrid->find_neighbors(A, neighbors);

		for (unsigned int i = 0; i < neighbors.size(); i++) {
			float t_collision = time_of_impact(*A, *neighbors[i], t_curr, m_time);

			if (t_collision < m_time) {
				m_time = t_collision; B = neighbors[i];
			}
		}
	} else {
		Point x = A->position();

		for (int i = -1; i <= 1; i++) {
			for (int j = -1; j <= 1; j++) {
				for (int k = -1; k <= 1; k++) {
					Vector shift(i*period[0], j*period[1], k*period[2]);

					neighbors.clear();
					A->set_position(x + shift);
					hgrid->find_neighbors(A, neighbors);

					for (unsigned int i = 0; i < neighbors.size(); i++) {
						float t_collision = time_of_impact(*A, *neighbors[i], t_curr, m_time);

						if (t_collision < m_time) {
							m_time = t_collision; B = neighbors[i]; m_shift = shift;
						}
					}
				}
			}
		}

		A->set_position(x); /* return particle to initial position */
	}

	if (B != NULL) { m_type = COLLISION; cB = B->collisions(); }
}

void Event::recompute() { predict(); }
