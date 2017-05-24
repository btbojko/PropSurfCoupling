#ifndef _PARTICLE_H
#define _PARTICLE_H

#include <vector>

#include "shape.h"
#include "state.h"
#include "transform.h"

class Particle : public State {
  public:
	Particle(Shape* shape, float density = 1.0)
		: m_shape(shape), m_growth_rate(1.0), m_collisions(0)
	{
		Vector zero(0.0, 0.0, 0.0);
		set_mass(density * volume());
		set_position(zero);
		set_velocity(zero);
		set_angular_velocity(zero);
		set_orientation(Quaternion::identity());
	}

	float mass()     const { return m_mass; }
	float inv_mass() const { return m_inv_mass; }

	void set_mass(float mass)
	{
		if (near_zero(mass)) {
			m_mass = 0.0; m_inv_mass = INFINITY;
		} else if (near_zero(1.0/mass)) {
			m_mass = INFINITY; m_inv_mass = 0.0;
		} else {
			m_mass = mass; m_inv_mass = 1.0/mass;
		}
	}

	Shape* shape() { return m_shape; }
	const Shape* shape() const { return m_shape; }
	void set_shape(Shape* s) { m_shape = s; }

	float growth_rate() const { return m_growth_rate; }
	void set_growth_rate(float rate) { m_growth_rate = rate; }

	float volume(float t = 1.0) const
	{
		return pow(t * m_growth_rate, 3) * m_shape->volume();
	}

	float bounding_radius(float t = 1.0) const
	{
		return t * m_growth_rate * m_shape->bounding_radius();
	}

	const Transform world_transform(float t) const {
		Transform T(orientation(t), position(t));
		T.scale(t * m_growth_rate);
		return T;
	}

	Matrix inverse_world_inertia()
	{
		Transform world(m_orientation, m_position);
		return m_inv_mass * world(m_shape->inverse_inertia());
	}

	Matrix inverse_world_inertia(float t)
	{
		Transform world(orientation(t), position(t));
		world.scale(t * m_growth_rate);
		return m_inv_mass * world(m_shape->inverse_inertia());
	}

	void apply_central_impulse(const Vector& impulse)
	{
		m_velocity += m_inv_mass * impulse;
	}

	void apply_angular_impulse(const Vector& torque)
	{
		m_ang_velocity += inverse_world_inertia() * torque;
	}

	void apply_impulse(const Vector& impulse, const Point& r)
	{
		if (m_inv_mass != 0.0) {
			apply_central_impulse(impulse);
			apply_angular_impulse(r ^ impulse);
		}
	}

	void unpack_state(float& inv_m, Point& p, Vector& v, Vector& w)
	{
		inv_m = m_inv_mass; p = m_position;
		v = m_velocity; w = m_ang_velocity;
	}

    Particle* prev() const { return m_prev; }
    Particle* next() const { return m_next; }
    void set_prev(Particle* p) { m_prev = p; }
    void set_next(Particle* p) { m_next = p; }

    unsigned int tag() const { return m_tag; }
    void set_tag(unsigned int id) { m_tag = id; }

    unsigned int get_event_id() const { return m_event_id; }
    void set_event_id(unsigned int id) { m_event_id = id; }

    unsigned long int hash() const { return m_hash; }
    void set_hash(unsigned long int hash) { m_hash = hash; }

    void collided() { m_collisions++; }
    unsigned long int collisions() const { return m_collisions; }

	void get_color(float& R, float& G, float& B, float& A)
	{
		R = r; G = g; B = b; A = a;
	}

	void set_color(float r, float g, float b, float a = 1.0)
	{
		assert(0.0 <= r && r <= 1.0);
		assert(0.0 <= g && g <= 1.0);
		assert(0.0 <= b && b <= 1.0);
		assert(0.0 <= a && a <= 1.0);
		m_color = _mm_set_ps(a, b, g, r);
	}

#ifdef HAVE_OPENGL
    void draw() const;
#endif

  private:
	Shape* m_shape;
	float m_growth_rate;
	float m_mass, m_inv_mass;
	unsigned int  m_tag;
	unsigned int  m_event_id;
	unsigned long int m_hash;
	unsigned long int m_collisions;
	Particle *m_prev, *m_next;
	union {
		__m128 m_color;
		struct { float r, g, b, a; };
	} __attribute__ ((aligned(16)));
};

extern std::vector<Particle*> particle;

#endif
