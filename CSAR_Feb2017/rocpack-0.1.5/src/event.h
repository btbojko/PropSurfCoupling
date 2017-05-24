#ifndef _EVENT_H
#define _EVENT_H

#include "vector.h"

class Particle;

enum EventType { INVALID, COLLISION, BOUNDARY_COLLISION, TRANSFER };

class Event {
  public:
	Event(unsigned int id, Particle* A);

	int type() const { return m_type; }
	float time() const { return m_time; }
	unsigned int id() const { return m_id; }
	unsigned int secondary_id() const;

	void print();
	void process();
	void predict();
	void recompute();

  private:
	float m_time;
	unsigned int m_id;
	Particle *A, *B;
	unsigned int cA, cB;
	Vector m_shift;
	EventType m_type;
};

#endif
