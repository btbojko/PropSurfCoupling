#ifndef _STATE_H
#define _STATE_H

#include "point.h"
#include "vector.h"
#include "quaternion.h"

class State {
  public:
	State() {}
	State(const Point p, const Vector v, const Vector w, const Quaternion q)
		: m_time(0.0), m_position(p), m_velocity(v),
		  m_ang_velocity(w), m_orientation(q) {}

	float time() const
	{
		return m_time;
	}

	const Point& position() const
	{
		return m_position;
	}

	const Point  position(float t) const
	{
		return m_position + (t - m_time) * m_velocity;
	}

	const Vector& velocity() const
	{
		return m_velocity;
	}

	const Vector& angular_velocity() const
	{
		return m_ang_velocity;
	}

	const Quaternion& orientation() const
	{
		return m_orientation;
	}

	const Quaternion orientation(float t) const
	{
		float s = norm(m_ang_velocity);
		Quaternion dq(m_ang_velocity/s, (t - m_time) * s);

		if (norm(dq) > EPSILON)
			return m_orientation ^ dq;
		else
			return m_orientation;
	}

	void set_time(float t)
	{
		m_time = t;
	}

	void set_position(const Point& p)
	{
		m_position = p;
	}

	void set_velocity(const Vector& v)
	{
		m_velocity = v;
	}

	void set_angular_velocity(const Vector& w)
	{
		m_ang_velocity = w;
	}

	void set_orientation(const Quaternion& q)
	{
		m_orientation = q;
		assert(abs(norm(q) - 1.0) < 1.0e-3);
	}

	void translate(const Vector& dx)
	{
		m_position += dx;
	}

	void rotate(const Quaternion& dq)
	{
		if (norm(dq) > EPSILON) {
			m_orientation ^= dq;
		}
		assert(abs(norm(m_orientation) - 1.0) < 1.0e-3);
	}

	void advance(float dt)
	{
		float s = norm(m_ang_velocity);
		m_time += dt;
		translate(dt * m_velocity);
		rotate(Quaternion(m_ang_velocity/s, dt*s));
	}

	void sync(float t)
	{
		float dt = t - m_time;
		float  s = norm(m_ang_velocity);
		m_time = t;
		translate(dt * m_velocity);
		rotate(Quaternion(m_ang_velocity/s, dt*s));
	}

  protected:
	float      m_time;
	Point      m_position;
	Vector     m_velocity;
	Vector     m_ang_velocity;
	Quaternion m_orientation;
};

#endif
