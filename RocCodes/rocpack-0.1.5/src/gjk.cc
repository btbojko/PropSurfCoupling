#include "convex.h"
#include "simplex.h"
#include "particle.h"
#include "transform.h"

static const float rel_error = 0.001;
static const float abs_error = 0.001;

bool contains(const Particle& A, const Point& x)
{

	Vector v = A.position() - x;

	/* early exit for when bounding volumes don't intersect */
	if (norm(v) > A.bounding_radius())
		return false;

	if (A.shape()->type() == SPHERE)
		return true;

	Convex& a = (Convex&)(*A.shape());
	Transform a2w = A.world_transform(1.0);
	Transform w2a = inverse(a2w);

	Simplex simplex; Point p, q = x; Vector w;

	do {
		p = a2w(a.support(w2a(-v)));
		w = p - q;

		if (v * w > 0.0)
			return false;

		if (simplex.contains(w))
			return false;

		simplex.add_vertex(w, p, q);

		if (!simplex.closest(v))
			return false;

	} while (!simplex.full() && norm(v) > EPSILON);

	return true;
}

bool intersect(const Particle& A, const Particle& B, float t)
{

	Vector v = A.position(t) - B.position(t);

	/* early exit for when bounding volumes don't intersect */
	if (norm(v) > A.bounding_radius(t) + B.bounding_radius(t))
		return false;

	int iterations = 0;
	Convex& a = (Convex&)(*A.shape());
	Convex& b = (Convex&)(*B.shape());
	Transform a2w = A.world_transform(t);
	Transform b2w = B.world_transform(t);
	Transform w2a = inverse(a2w), w2b = inverse(b2w);

	Simplex simplex; Point p, q; Vector w;

	do {
		p = a2w(a.support(w2a(-v)));
		q = b2w(b.support(w2b( v)));
		w = p - q;

		if (v * w > 0.0)
			return false;

		if (simplex.contains(w))
			return false;

		simplex.add_vertex(w, p, q);

		if (!simplex.closest(v))
			return false;

	} while (!simplex.full() &&
	         iterations++ < 50 && norm(v) > abs_error);

	return true;
}

float distance(const Particle& A, const Particle& B, float t)
{
	Convex& a = (Convex&)(*A.shape());
	Convex& b = (Convex&)(*B.shape());
	Transform a2w = A.world_transform(t);
	Transform b2w = B.world_transform(t);
	Transform w2a = inverse(a2w), w2b = inverse(b2w);

	Simplex simplex; Point p, q; Vector v, w;

	v = A.position(t) - B.position(t);

	float dist = norm(v), mu = 0.0;

	do {
		p = a2w(a.support(w2a(-v)));
		q = b2w(b.support(w2b( v)));
		w = p - q;

		set_max(mu, (v * w) / dist);

		if (simplex.contains(w))
			return dist;

		simplex.add_vertex(w, p, q);

		if (dist - mu < rel_error * dist)
			return dist;

		if (!simplex.closest(v))
			return dist;

		dist = norm(v);

	} while (!simplex.full() && dist > abs_error);

	return 0.0;
}

bool closest_points(const Particle& A, const Particle& B, Point& pA, Point& pB, float t)
{
	Convex& a = (Convex&)(*A.shape());
	Convex& b = (Convex&)(*B.shape());
	Transform a2w = A.world_transform(t);
	Transform b2w = B.world_transform(t);
	Transform w2a = inverse(a2w), w2b = inverse(b2w);

	Simplex simplex; Point p, q; Vector v, w;

	v = A.position(t) - B.position(t);

	float dist = norm(v), mu = 0.0;

	do {
		p = a2w(a.support(w2a(-v)));
		q = b2w(b.support(w2b( v)));
		w = p - q;

		set_max(mu, (v * w) / dist);

		if (simplex.contains(w))
			break;

		simplex.add_vertex(w, p, q);

		if (dist - mu < rel_error * dist)
			break;

		if (!simplex.closest(v))
			break;

		dist = norm(v);

	} while (!simplex.full() && dist > abs_error);

	assert(!simplex.empty());

	simplex.compute_points(pA, pB);

	return !simplex.full();
}

bool contact(const Particle& A, const Particle& B, Point& pA, Point& pB, Vector& N, float t)
{
	closest_points(A, B, pA, pB, t);

	N = (pA - pB).normalized();

	/* catch NaNs */
	assert(norm(pA) > 0.0);
	assert(norm(pB) > 0.0);
	assert(norm( N) > 0.0);

	return true;
}
