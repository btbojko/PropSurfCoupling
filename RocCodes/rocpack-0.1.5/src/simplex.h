#ifndef _SIMPLEX_H
#define _SIMPLEX_H

#include "point.h"
#include "vector.h"

struct closest_result {
	Point closest_point;
	float barycentric_coords[4];
	bool  degenerate, valid;

	void reset()
	{
		set_barycentric_coordinates();
		valid = degenerate = false;
	}

	void set_barycentric_coordinates(float a = 0.0, float b = 0.0,
	                                 float c = 0.0, float d = 0.0)
	{
		barycentric_coords[0] = a;
		barycentric_coords[1] = b;
		barycentric_coords[2] = c;
		barycentric_coords[3] = d;

		valid = (a >= 0.0 && b >= 0.0 && c >= 0.0 && d >= 0.0);
	}
};

class Simplex {
  public:
	Simplex() { reset(); }

	void reset();

	bool full()  const { return n_vertices == 4; }
	bool empty() const { return n_vertices == 0; }

	int  vertices() const { return n_vertices; }

	void add_vertex(const Vector& w, const Point& p, const Point& q);
	void remove_vertex(int index);
	bool contains(const Vector& w);

	void reduce();
	bool update();
	bool degenerate() const;

	bool closest(Vector& v);
	void backup_closest(Vector& v);

	void compute_points(Point& p, Point& q);

  private:
	bool out_of_date;

	int n_vertices;

	Point  m_p[4]; /* points of object A in world coordinates */
	Point  m_q[4]; /* points of object B in world coordinates */
	Vector m_w[4]; /* points of  (A - B) in world coordinates */

	Point  m_cached_p;
	Point  m_cached_q;
	Vector m_cached_v;
	Vector m_last_w;

	struct closest_result m_closest;
};

#endif
