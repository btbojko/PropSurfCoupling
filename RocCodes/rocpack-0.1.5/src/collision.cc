#include <cstdio>
#include <algorithm>

#include "stime.h"
#include "gjk.h"
#include "queue.h"
#include "hgrid.h"
#include "particle.h"
#include "boundary.h"
#include "collision.h"

extern float target_fraction;
static const float collision_margin = 5.0e-3;
static const float precision_limit  = 1.0e-3;

void contact_network(Particle& p, FILE* out)
{
	std::vector<float> dist;
	std::vector<Particle*> neighbors;

	if (boundary->is_periodic()) {
		Point x = p.position();
		for (int i = -1; i <= 1; i++) {
			for (int j = -1; j <= 1; j++) {
				for (int k = -1; k <= 1; k++) {
					Vector shift(i*period[0], j*period[1], k*period[2]);

					neighbors.clear();
					p.set_position(x + shift);
					hgrid->find_neighbors(&p, neighbors);

					for (unsigned int i = 0; i < neighbors.size(); i++) {
						if (p.shape()->type() == SPHERE &&
							neighbors[i]->shape()->type() == SPHERE) {
							float d = (p.position() - neighbors[i]->position()).norm() -
								t_curr * (p.growth_rate() + neighbors[i]->growth_rate());
							dist.push_back(d);
						} else {
							dist.push_back(distance(p, *neighbors[i], t_curr));
						}
					}
				}
			}
		}

		p.set_position(x); /* return particle to initial position */

	} else {
		neighbors.clear();
		hgrid->find_neighbors(&p, neighbors);

		for (unsigned int i = 0; i < neighbors.size(); i++) {
			dist.push_back(distance(p, *neighbors[i], t_curr));
		}
	}

	std::sort(dist.begin(), dist.end());

	for (unsigned int i = 0; i < dist.size(); i++) {
		fprintf(out, " %.3e", dist[i] / p.bounding_radius(t_curr));
	}

	fprintf(out, "\n");
}

float sphere_time_of_impact(const Particle& A, const Particle& B)
{
	float a, b, c, d, xv, x2, v2, r, r2, vrel;
	Vector x, v;

	x  = B.position(t_curr) - A.position(t_curr);
	v  = B.velocity() - A.velocity();
	r  = A.bounding_radius() + B.bounding_radius();
	xv = x*v; x2 = x*x; v2 = v*v; r2 = r*r; vrel = xv/sqrt(x2) - r;

	c  = x2 - t_curr*t_curr*r2;

	if (c > 0.0) { /* spheres are not touching */
		/* solve quadratic equation for the minimum distance */
		a = v2 - r2; b = xv - r2*t_curr; d = b*b - a*c;

        if (vrel >= -0.0 || d < 0.0)
			return INFINITY;

		float q, t;

		if (b >= 0.0) {
			q = -(b + sqrt(d)); t = t_curr + q/a;
		} else {
			q = -(b - sqrt(d)); t = t_curr + c/q;
		}

		/* intersection happened only in the past */
		if (t < t_curr) return INFINITY;

		return t;

	} else { /* spheres are touching */
		if (vrel <= 0.0) return t_curr; else return INFINITY;
    }
}

void sphere_collision_response(Particle& A, Particle& B)
{
	float e = 1.0; /* restitution coefficient */
	float inv_mA, inv_mB;
	Point  xA, xB, pA, pB;
	Vector vA, wA, vB, wB;

	A.unpack_state(inv_mA, xA, vA, wA);
	B.unpack_state(inv_mB, xB, vB, wB);
	float g = A.bounding_radius() + B.bounding_radius();

	Vector  N = (xA - xB).normalized(); /* points from B to A */
	Vector  v = vA - vB;

	float vrel = -N*v;

	vrel = vrel <= 0.0 ? (1+EPSILON) * g : vrel + g;

	Vector impulse = ((1.0 + e) * (vrel / (inv_mA + inv_mB))) * N;

	A.set_velocity(vA + inv_mA * impulse);
	B.set_velocity(vB - inv_mB * impulse);
}

void overlap_interval(const Particle& A,
                      const Particle& B, float& tmin, float& tmax)
{
	float a, b, c, d, q, g, g0;

	Vector x = B.position(t_curr) - A.position(t_curr);
	Vector v = B.velocity() - A.velocity();
	g  = A.bounding_radius() + B.bounding_radius();
	g0 = g*t_curr + 1.2*collision_margin;
	a = v*v - g*g; b = x*v - g*g0; c = x*x - g0*g0;
	d = b*b - a*c;

	/* no intersection */
	if (d < 0.0) {
		tmin = tmax = INFINITY;
		return;
	}

	/* some intersection */
	if (b >= 0.0) {
		q = -(b + sqrt(d));
		tmin = t_curr + q/a;
		tmax = t_curr + c/q;
	} else {
		q = -(b - sqrt(d));
		tmin = t_curr + c/q;
		tmax = t_curr + q/a;
	}

	/* growing faster than separating */
	if (tmax < tmin) tmax = INFINITY;
}

float relative_velocity(const Particle& A, const Particle& B, float t)
{
	Point  xA = A.position(t);
	Point  xB = B.position(t);
	Vector vA = A.velocity();
	Vector vB = B.velocity();
	Vector wA = A.angular_velocity();
	Vector wB = B.angular_velocity();

	Point pA, pB; Vector N;

	float vrel, dist = distance(A, B, t);

	if (dist > precision_limit) {
		contact(A, B, pA, pB, N, t);
		Vector rA = pA - xA;
		Vector rB = pB - xB;

		if (N * rB < 0.0) N = -N;

		vrel = N * ((vA + (wA^rA)) - (vB + (wB^rB)));
	} else {
		vrel = (xA - xB) * (vA - vB) / norm(xA - xB);
	}

	return vrel;
}

bool collision_check(const Particle& A, const Particle& B, float t)
{
	float dist1 = distance(A, B, t);
	float dist2 = distance(A, B, t + 0.00001);

	return (dist2 < dist1);
}

float bisect_time_of_impact(const Particle& A, const Particle& B, float tmin, float tmax)
{

	if (distance(A, B, tmin) < collision_margin)
		return tmin;

	while (distance(A, B, tmin) > 1.01 * collision_margin &&
	       !nearly_equal(tmin, tmax)) {
		float tmid = 0.5 * (tmin + tmax);

		if (distance(A, B, tmid) > collision_margin)
			tmin = tmid;
		else
			tmax = tmid;
	}

	return tmin;
}

float search_time_of_impact(const Particle& A, const Particle& B,
                            float tmin, float tmax, unsigned int depth = 0)
{
	static const unsigned int max_depth = 8;

	if (distance(A, B, tmin) < collision_margin)
		return depth == 0 ? tmin : INFINITY;

	float t = 0.5 * (tmin + tmax);

	if (distance(A, B, t) < collision_margin) {
		float t1 = bisect_time_of_impact(A, B, tmin, t);
		float t2 = search_time_of_impact(A, B, tmin, t1, depth+1);
		return min<float>(t1, t2);
	}

	if (nearly_equal(tmin, tmax) || depth >= max_depth)
		return INFINITY;

	float t1 = search_time_of_impact(A, B, tmin, t, depth+1);

	if (t1 < INFINITY)
		return t1;
	else
		return search_time_of_impact(A, B, t, tmax, depth+1);
}

float convex_time_of_impact(const Particle& A, const Particle& B, float tmin, float tmax)
{
	float t1, t2;

	overlap_interval(A, B, t1, t2);

	if (t1 > t_stop || t1 > tmax || t2 < tmin)
		return INFINITY;

	set_max(t1, tmin); set_min(t2, tmax);

	if (distance(A, B, t1) < 1.001 * collision_margin) {
		if (collision_check(A, B, t1)) {
			return t1;
		} else {
			return convex_time_of_impact(A, B, t1 + EPSILON, t2);
		}
	}

	return search_time_of_impact(A, B, t1, t2, 0);
}

void convex_collision_response(Particle& A, Particle& B)
{
	static const float e = 1.0;
	float  inv_mA, inv_mB;
	float  gA = A.growth_rate();
	float  gB = B.growth_rate();
	Point  xA, xB, pA, pB;
	Vector vA, wA, vB, wB, N;

	A.unpack_state(inv_mA, xA, vA, wA);
	B.unpack_state(inv_mB, xB, vB, wB);

	contact(A, B, pA, pB, N, t_curr);

	assert(N * (xA - xB) > 0.0);

	Point p = (pA + pB) / 2.0;
	Vector rA = p - xA, kA = rA ^ N;
	Vector rB = p - xB, kB = rB ^ N;
	Vector uA = A.inverse_world_inertia() * kA;
	Vector uB = B.inverse_world_inertia() * kB;

	float vrel  = N * ((vA+(wA^rA)) - (vB+(wB^rB)));
	float denom = (inv_mA + inv_mB + kA*uA + kB*uB);

	if (vrel >= -(gA + gB))
		vrel = (vrel >= 0.0) ? -(1.01 * (gA + gB)):
							    (vrel - (gA + gB));

	float f = -(1.0 + e) * vrel / denom;

	Vector impulse = f * N;

	A.apply_impulse( impulse, rA);
	B.apply_impulse(-impulse, rB);
}

#ifndef EDMC

float time_of_impact(const Particle& A, const Particle& B, float tmin, float tmax)
{
	if (A.shape()->type() == SPHERE && B.shape()->type() == SPHERE) {
		return sphere_time_of_impact(A, B);
	} else {
		return convex_time_of_impact(A, B, tmin, tmax);
	}
}

void collision_response(Particle& A, Particle& B)
{
	if (A.shape()->type() == SPHERE && B.shape()->type() == SPHERE) {
		sphere_collision_response(A, B);
	} else {
		convex_collision_response(A, B);
	}
}

#else

bool check_overlap(Particle& p)
{
	static std::vector<Particle*> neighbors;

	if (boundary->is_periodic()) {
		Point x = p.position();
		for (int i = -1; i <= 1; i++) {
			for (int j = -1; j <= 1; j++) {
				for (int k = -1; k <= 1; k++) {
					Vector shift(i*period[0], j*period[1], k*period[2]);

					neighbors.clear();
					p.set_position(x + shift);
					hgrid->find_neighbors(&p, neighbors);

					for (unsigned int i = 0; i < neighbors.size(); i++) {
						if (intersect(p, *neighbors[i], t_curr)) {
							p.set_position(x);
							return false;
						}
					}
				}
			}
		}

		p.set_position(x); /* return particle to initial position */

	} else {
		neighbors.clear();
		hgrid->find_neighbors(&p, neighbors);

		for (unsigned int i = 0; i < neighbors.size(); i++) {
			if (intersect(p, *neighbors[i], t_curr))
				return false;
		}
	}

	return true;
}

bool move_particle(Particle& p)
{
	int n = 1, tries;
	Point x = p.position();
	Quaternion q = p.orientation();

	float f = target_fraction * pow(t_curr/t_stop, 3.0);

	tries = (int) (100.0 / pow(1.1 - f, 2.0));
	float coeff = pow(1.0 - f, 3.0) + EPSILON;

#ifndef LATTICE_PACKING
	if (p.shape()->type() != SPHERE) {
		do {
			float s = coeff / n;
			Vector d = Vector::random();
			Quaternion dq = Quaternion(d, s);
			p.set_orientation(q ^ dq);
		} while (!check_overlap(p) && ++n < tries);

		n = 1;
	}
#endif

	do {
		float s = coeff / n;
		Vector dx = x + s * Vector::random();
		p.set_position(dx);

		/* keep particle inside boundary */
		if (boundary->is_periodic())
			((PeriodicBoundary*)boundary)->reposition(p);

		hgrid->rehash(&p);
	} while (!check_overlap(p) && ++n < tries);

	/* revert changes if moving failed */

	if (n >= tries) {
		p.set_position(x);
		p.set_orientation(q);
		hgrid->rehash(&p);
		return false;
	}

	return true;
}

float time_of_impact(const Particle& A,
                     const Particle& B, float tmin, float tmax)
{
	unsigned int n = 0;

	tmin = t_curr;
	tmax = t_stop;

	if (intersect(A, B, t_curr)) /* because of floating point errors */
		return t_curr;

	if(!intersect(A, B, tmax))
		return INFINITY;

	while (tmax - tmin > EPSILON && n++ < 15) {
		float tmid = (tmax + tmin)/2.0;
		if (intersect(A, B, tmid))
			tmax = tmid;
		else
			tmin = tmid;
	}

	return tmin;
}

void collision_response(Particle& A, Particle& B)
{
	static std::vector<Particle*> neighbors;

	if(!move_particle(A) && !move_particle(B)) {
			neighbors.clear();
			hgrid->find_neighbors(&A, neighbors);
			for (unsigned int i = 0; i < neighbors.size(); i++) {
				move_particle(*neighbors[i]);
				event_queue->update_id(neighbors[i]->get_event_id());
			}

			neighbors.clear();
			hgrid->find_neighbors(&B, neighbors);
			for (unsigned int i = 0; i < neighbors.size(); i++) {
				move_particle(*neighbors[i]);
				event_queue->update_id(neighbors[i]->get_event_id());
			}

			move_particle(A);
			move_particle(B);
	}
}
#endif
