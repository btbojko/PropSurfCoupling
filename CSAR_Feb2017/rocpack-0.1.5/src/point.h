#ifndef _POINT_H
#define _POINT_H

#include "vector.h"

class Point : public Vector {
  public:
	Point() { data = _mm_setzero_ps(); }

	Point(__m128 m) : Vector(m) {}

	Point(float x, float y, float z) : Vector(x, y, z) {}

	Point(const Vector& v) : Vector(v) {}
};

inline float min(const Point& p) { return p.min(); }
inline float max(const Point& p) { return p.max(); }
inline Point abs(const Point& p) { return p.abs(); }
inline Point sgn(const Point& p) { return p.sgn(); }

inline float norm(const Point& p) { return p.norm(); }
inline float rnorm(const Point& p) { return p.rnorm(); }
inline float norm2(const Point& p) { return p.norm2(); }

inline Point operator+(const Point& p, const Vector& v)
{
	return Point(_mm_add_ps(p.get_value(), v.get_value()));
}

inline Point operator-(const Point& p, const Vector& v)
{
	return Point(_mm_sub_ps(p.get_value(), v.get_value()));
}

inline Vector operator-(const Point& p1, const Point& p2)
{
	return Vector(_mm_sub_ps(p1.get_value(), p2.get_value()));
}

inline Point min(const Point& p1, const Point& p2)
{
	return Point(_mm_min_ps(p1.get_value(), p2.get_value()));
}

inline Point max(const Point& p1, const Point& p2)
{
	return Point(_mm_max_ps(p1.get_value(), p2.get_value()));
}

inline float distance(const Point& a, const Point& b)
{
	return norm(a - b);
}

inline float distance2(const Point& a, const Point& b)
{
	return norm2(a - b);
}

#endif
