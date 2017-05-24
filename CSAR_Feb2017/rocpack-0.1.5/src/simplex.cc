#include "scalar.h"
#include "simplex.h"

void Simplex::reset()
{
	out_of_date = true; n_vertices = 0; m_closest.reset();
	m_last_w.set_value(INFINITY, INFINITY, INFINITY);
}

void Simplex::add_vertex(const Vector& w, const Point& p, const Point& q)
{
	assert(!full());

	out_of_date = true;

	m_last_w = w;

	m_p[n_vertices] = p;
	m_q[n_vertices] = q;
	m_w[n_vertices] = w;

	n_vertices++;
}

void Simplex::remove_vertex(int index)
{
	assert(n_vertices > 0);

	n_vertices--;
	m_p[index] = m_p[n_vertices];
	m_q[index] = m_q[n_vertices];
	m_w[index] = m_w[n_vertices];
}

bool Simplex::contains(const Vector& w)
{
	const float threshold = EPSILON;

	for (int i = 0; i < n_vertices; i++)
		if (norm(m_w[i] - w) < threshold)
			return true;

	if (norm(m_last_w - w) < EPSILON)
		return true;

	return false;
}

void Simplex::reduce()
{
	for (int i = 3; i >= 0; --i) {
		if (n_vertices > i && m_closest.barycentric_coords[i] == 0.0)
			remove_vertex(i);
	}
}

bool Simplex::closest(Vector& v)
{
	bool status = update();
	if (m_closest.valid)
		v = m_cached_v;
	return status;
}

void Simplex::backup_closest(Vector& v)
{
	v = m_cached_v;
}

void Simplex::compute_points(Point& p, Point& q)
{
	update(); p = m_cached_p; q = m_cached_q;
}

void closest_point_triangle(const Point& a, const Point& b, const Point& c,
                            struct closest_result& result);

bool closest_point_tetrahedron(const Point& a, const Point& b,
                               const Point& c, const Point& d,
							   struct closest_result&  result);

bool Simplex::update()
{
	if (out_of_date) {

		switch (n_vertices) {

			case 0: /* simplex is empty */
				m_closest.reset();
				break;

			case 1: /* simplex is a single point */
				m_closest.reset();
				m_closest.closest_point = m_w[0];
				m_closest.set_barycentric_coordinates(1.0, 0.0, 0.0, 0.0);
				break;

			case 2: /* simplex is a line segment */
			{
				Vector ab = m_w[1] - m_w[0];
				float t = -(m_w[0] * ab) / norm2(ab);

				clamp(t, 0.0f, 1.0f);

				m_closest.closest_point = m_w[0] + t * ab;
				m_closest.set_barycentric_coordinates(1-t, t, 0.0, 0.0);
				break;
			}

			case 3: /* simplex is a triangle */
				closest_point_triangle(m_w[0], m_w[1], m_w[2], m_closest);
				break;

			case 4: /* simplex is a tetrahedron */
			{
				bool separated =
					closest_point_tetrahedron(m_w[0], m_w[1],
				                              m_w[2], m_w[3], m_closest);

				if (!separated) {
					if (!m_closest.degenerate) {
						m_closest.valid = true;
						m_cached_v.set_value(0.0, 0.0, 0.0);
					}

					/* do not update cached closest points if simplex is
					 * degenerate, use backup instead */
					return m_closest.valid;
				}

				break;
			}

			default:
				m_closest.valid = false;
		}

		m_cached_p.set_value(0.0, 0.0, 0.0);
		m_cached_q.set_value(0.0, 0.0, 0.0);

		for (int i = 0; i < n_vertices; i++) {
			m_cached_p += m_closest.barycentric_coords[i] * m_p[i];
			m_cached_q += m_closest.barycentric_coords[i] * m_q[i];
		}

		m_cached_v = m_cached_p - m_cached_q;

		reduce();

		out_of_date = false;
	}

	return m_closest.valid;
}

void closest_point_triangle(const Point& a, const Point& b, const Point& c,
                            struct closest_result& result)
{
	Vector ab = b - a, ac = c - a;
	float d1 = -a * ab, d2 = -a * ac;

	result.reset();

	/* closest feature is vertex a */

	if (d1 <= 0.0 && d2 <= 0.0) {
		result.closest_point = a;
		result.set_barycentric_coordinates(1.0, 0.0, 0.0, 0.0);
		return;
	}

	float d3 = -b * ab, d4 = -b * ac;

	/* closest feature is vertex b */

	if (d3 >= 0.0 && d4 <= d3) {
		result.closest_point = b;
		result.set_barycentric_coordinates(0.0, 1.0, 0.0, 0.0);
		return;
	}

	float vc = d1*d4 - d3*d2;

	/* closest feature is edge ab */
	if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0) {
		float t = (d1 / (d1 - d3));
		result.closest_point = a + t * ab;
		result.set_barycentric_coordinates(1-t, t, 0.0, 0.0);
		return;
	}

	float d5 = -c * ab, d6 = -c * ac;

	/* closest feature is vertex c */

	if (d6 >= 0.0 && d5 <= d6) {
		result.closest_point = c;
		result.set_barycentric_coordinates(0.0, 0.0, 1.0, 0.0);
		return;
	}

	float vb = d5*d2 - d1*d6;

	/* closest feature is edge ac */

	if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0) {
		float t = (d2 / (d2 - d6));
		result.closest_point = a + t * ac;
		result.set_barycentric_coordinates(1-t, 0.0, t, 0.0);
		return;
	}

	float va = d3*d6 - d5*d4;

	/* closest feature is edge bc */

	if (va <= 0.0 && (d4-d3) >= 0.0 && (d5-d6) >= 0.0) {
		float t = ((d4 - d3) / ((d4 - d3) + (d5 - d6)));
		result.closest_point = b + t * (c - b);
		result.set_barycentric_coordinates(0.0, 1-t, t, 0.0);
		return;
	}

	/* closest point lies on the triangle */

	float d = 1.0 / (va + vb + vc), v = vb * d, w = vc * d;
	result.closest_point = a + v*ab + w*ac;
	result.set_barycentric_coordinates(1-v-w, v, w, 0.0);
	return;
}

inline int out_of_plane(const Vector& a, const Vector& b, const Vector& c, const Vector& d)
{
	Vector normal = (b-a) ^ (c-a);

	float  sign1 = ( -a) * normal;
	float  sign2 = (d-a) * normal;

	if (sign2 * sign2 < 5.0e-14)
		return -1;

	return sign1 * sign2 < 0.0;
}

/* returns false if there is an intersection */

bool closest_point_tetrahedron(const Point& a, const Point& b,
                               const Point& c, const Point& d,
							   struct closest_result&  result)
{
	result.reset();
	result.closest_point.set_value(0.0, 0.0, 0.0);

	int out_abc = out_of_plane(a, b, c, d);
	int out_acd = out_of_plane(a, c, d, b);
	int out_adb = out_of_plane(a, d, b, c);
	int out_bcd = out_of_plane(b, d, c, a);

	/* points are affine dependent / degenerate */

	if (out_abc < 0 || out_acd < 0 || out_adb < 0 || out_bcd < 0) {
		result.degenerate = true;
		return false;
	}

	/* origin is inside the tetrahedron, i.e., A and B intersect */

	if (!out_abc && !out_acd && !out_adb && !out_bcd)
		return false;

	/* no intersection, loop over faces to find closest point */

	struct closest_result tmp_result;
	float u, v, w, dist, closest_dist = INFINITY;

	if (out_abc) {
		closest_point_triangle(a, b, c, tmp_result);

		dist = norm(tmp_result.closest_point);

		if (dist < closest_dist) {
			closest_dist = dist;
			u = tmp_result.barycentric_coords[0];
			v = tmp_result.barycentric_coords[1];
			w = tmp_result.barycentric_coords[2];
			result.closest_point = tmp_result.closest_point;
			result.set_barycentric_coordinates(u, v, w, 0.0);
		}
	}

	if (out_acd) {
		closest_point_triangle(a, c, d, tmp_result);

		dist = norm(tmp_result.closest_point);

		if (dist < closest_dist) {
			closest_dist = dist;
			u = tmp_result.barycentric_coords[0];
			v = tmp_result.barycentric_coords[1];
			w = tmp_result.barycentric_coords[2];
			result.closest_point = tmp_result.closest_point;
			result.set_barycentric_coordinates(u, 0.0, v, w);
		}
	}

	if (out_adb) {
		closest_point_triangle(a, d, b, tmp_result);

		dist = norm(tmp_result.closest_point);

		if (dist < closest_dist) {
			closest_dist = dist;
			u = tmp_result.barycentric_coords[0];
			v = tmp_result.barycentric_coords[1];
			w = tmp_result.barycentric_coords[2];
			result.closest_point = tmp_result.closest_point;
			result.set_barycentric_coordinates(u, w, 0.0, v);
		}
	}

	if (out_bcd) {
		closest_point_triangle(b, c, d, tmp_result);

		dist = norm(tmp_result.closest_point);

		if (dist < closest_dist) {
			closest_dist = dist;
			u = tmp_result.barycentric_coords[0];
			v = tmp_result.barycentric_coords[1];
			w = tmp_result.barycentric_coords[2];
			result.closest_point = tmp_result.closest_point;
			result.set_barycentric_coordinates(0.0, u, v, w);
		}
	}

	return true;
}
