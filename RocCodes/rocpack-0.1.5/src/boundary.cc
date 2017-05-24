#include "stime.h"
#include "hgrid.h"
#include "boundary.h"
#include "particle.h"

#ifdef HAVE_OPENGL

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#endif

Vector period;
Boundary* boundary;

Point Boundary::get_random_position() const
{
	Point p;

	do {
		float rngx = 0.5 - drand48();
		float rngy = 0.5 - drand48();
		float rngz = 0.5 - drand48();
		Vector rng(rngx, rngy, rngz);
		p = elem_mul(rng, dimensions());
	} while ((*this)(p) >= -10 * EPSILON);

	return p;
}

float Boundary::predict_collision(const Particle& p) const
{
	float dt_min = 0.0, dt_max = t_stop;
	float distance_to_boundary =
		-((*this)(p.position(t_curr))
		+ p.bounding_radius(t_curr));

	assert(distance_to_boundary > 0.0);

	while (dt_max - dt_min > EPSILON) {
		float dt = 0.5 * (dt_max + dt_min);

		distance_to_boundary =
			-((*this)(p.position(t_curr + dt))
			+ p.bounding_radius(t_curr + dt));

		if (distance_to_boundary > 0.001) {
			dt_min = dt;
		} else {
			dt_max = dt;
		}
	}

	return t_curr + dt_min;
}

void Boundary::collision_response(Particle& p) const
{
	Point  x = p.position();
	Vector v = p.velocity();
	float  g = p.growth_rate();
	float vrel = v * gradient(x) + g;
	float deltav = (vrel >= g ? -2.0001 * vrel : -1.0001 * g);

	p.apply_central_impulse(p.mass() * deltav * gradient(x));
}

float PeriodicBoundary::operator()(const Vector& x) const
{
	return max(abs(x) - size/2.0);
}

Vector PeriodicBoundary::gradient(const Vector& x) const
{
	int max = abs(elem_div(x, size)).max_element();
	Vector v(0.0, 0.0, 0.0); v[max] = sgn(x[max]);
	return v;
}

float PeriodicBoundary::volume() const
{
	return size[0]*size[1]*size[2];
}

#ifdef HAVE_OPENGL
void PeriodicBoundary::draw() const
{
	glPushMatrix();
	glScalef(size[0], size[1], size[2]);
	glutWireCube(1.0);
	glPopMatrix();
}
#endif

float PeriodicBoundary::predict_collision(const Particle& p) const
{
	float dt[3] = {t_stop, };
	Vector x = p.position();
	Vector v = p.velocity();

	for (int i = 0; i < 3; i++) {
		if (v[i] >= 0.0) {
			assert(size[i]/2.0 - x[i] >= 0.0);
			dt[i] = abs(( size[i]/2.0 + 10*EPSILON - x[i]) / v[i]);
		} else {
			assert(-size[i]/2.0 - x[i] <= 0.0);
			dt[i] = abs((-size[i]/2.0 - 10*EPSILON - x[i]) / v[i]);
		}
	}

	return t_curr + max<float>(min(dt[0], dt[1], dt[2]), EPSILON);
}

void PeriodicBoundary::collision_response(Particle& p) const
{
	Point  x = p.position();
	Vector v = p.velocity();

	for (int i = 0; i < 3; i++) {
		if (v[i] < 0.0 && x[i] < -size[i]/2.0)
			x[i] += size[i];
		else if (v[i] > 0.0 && x[i] > size[i]/2.0)
			x[i] -= size[i];
	}

	p.set_position(x);
	hgrid->rehash(&p);
}

/* keep particle within the periodic cell */
void PeriodicBoundary::reposition(Particle& p) const
{
	Point x = p.position();

	while (abs(x[0]) > size[0]/2.0 ||
	       abs(x[1]) > size[1]/2.0 ||
		   abs(x[2]) > size[2]/2.0)
	{
		for (int i = 0; i < 3; i++) {
			if (x[i] < -size[i]/2.0) {
				x[i] += size[i];
			} else if (x[i] > size[i]/2.0) {
				x[i] -= size[i];
			}
		}
	}

	p.set_position(x);
	hgrid->rehash(&p);
}

float CubicBoundary::operator()(const Vector& x) const
{
	return max(abs(x)) - a/2.0;
}

Vector CubicBoundary::gradient(const Vector& x) const
{
	int max = abs(x).max_element();
	Vector v(0.0, 0.0, 0.0);
	v[max] = sgn(x[max]);
	return v;
}

float CubicBoundary::predict_collision(const Particle& p) const
{
	float dt = t_stop;
	Vector x = p.position();
	Vector v = p.velocity();
	float  g = p.growth_rate();
	float tolerance = 0.0005;

	for (int i = 0; i < 3; i++) {
		if (v[i] + g > 0.0) {
			float delta = a/2.0 - x[i] - g*t_curr - 3*tolerance;
			set_min<float>(dt, delta / (v[i] + g));
		}

		if (v[i] - g < 0.0) {
			float delta = g*t_curr - a/2.0 - x[i] + 3*tolerance;
			set_min<float>(dt, delta / (v[i] - g));
		}
	}

	return t_curr + max<float>(dt, 0.0);
}

void CubicBoundary::collision_response(Particle& p) const
{
	Vector x = p.position();
	Vector v = p.velocity();
	float  g = p.growth_rate();

	for (int i = 0; i < 3; i++) {
		if (a/2.0 - abs(x[i]) - g*t_curr < 0.005) {
			if (x[i] >= 0.0 && v[i] + g >= 0.0) {
				v[i] = 0.9*v[i] <  g ? -1.001 * g : -0.900 * v[i];
			} else if (x[i] <= 0.0 && v[i] - g <= 0.0) {
				v[i] = 0.9*v[i] > -g ?  1.001 * g : -0.900 * v[i];
			}
		}
	}

	p.set_velocity(v);
}

float CubicBoundary::volume() const
{
	return a * a * a;
}

#ifdef HAVE_OPENGL
void CubicBoundary::draw() const
{
	glPushMatrix();
	glutWireCube(a); /* a is half the edge size */
	glPopMatrix();
}
#endif


float BoxBoundary::operator()(const Vector& x) const
{
	return max(abs(x) - size/2.0);
}

Vector BoxBoundary::gradient(const Vector& x) const
{
	int max = abs(elem_div(x, size)).max_element();
	Vector v(0.0, 0.0, 0.0);
	v[max] = sgn(x[max]);
	return v;
}

float BoxBoundary::volume() const
{
	return size[0] * size[1] * size[2];
}

#ifdef HAVE_OPENGL
void BoxBoundary::draw() const
{
	glPushMatrix();
	glScaled(size[0], size[1], size[2]);
	glutWireCube(1.0);
	glPopMatrix();
}
#endif

Point BoxBoundary::get_random_position() const
{
	return (1.0 - EPSILON) * Vector(size[0] * (0.5-drand48()),
	                                size[1] * (0.5-drand48()),
									size[2] * (0.5-drand48()));
}

float BoxBoundary::predict_collision(const Particle& p) const
{
	float dt = t_stop;
	Vector x = p.position();
	Vector v = p.velocity();
	float  g = p.growth_rate();
	float tolerance = 0.001;

	for (int i = 0; i < 3; i++) {
		if (v[i] + g >= 0.0)
			set_min<float>(dt, (size[i]/2.0 - x[i] - g*t_curr - tolerance) / (v[i] + g));

		if (v[i] - g <= 0.0)
			set_min<float>(dt, (g*t_curr - size[i]/2.0 - x[i] + tolerance) / (v[i] - g));
	}

	return t_curr + max<float>(dt, 0.0);
}

void BoxBoundary::collision_response(Particle& p) const
{
	Vector x = p.position();
	Vector v = p.velocity();
	float  g = p.growth_rate();

	for (int i = 0; i < 3; i++) {
		if (size[i]/2.0 - abs(x[i]) - g*t_curr < 0.005) {
			if (x[i] >= 0.0 && v[i] + g >= 0.0) {
				v[i] = v[i] <  g ? -1.001 * g : -v[i];
			} else if (x[i] <= 0.0 && v[i] - g <= 0.0) {
				v[i] = v[i] > -g ?  1.001 * g : -v[i];
			}
		}
	}

	p.set_velocity(v);
}


float SphericBoundary::operator()(const Vector& x) const
{
	return norm(x) - r;
}

Vector SphericBoundary::gradient(const Vector& x) const
{
	return x / norm(x);
}

float SphericBoundary::volume() const
{
	return 4.0 / 3.0 * M_PI * r * r * r;
}

#ifdef HAVE_OPENGL
void SphericBoundary::draw() const
{
	glPushMatrix();
	glutWireSphere(r, 72, 72);
	glPopMatrix();
}
#endif


float CylindricBoundary::operator()(const Vector& x) const
{
	return max<float>(sqrt(x[0]*x[0] + x[2]*x[2]) - r, abs(x[1]) - h/2.0);
}

Vector CylindricBoundary::gradient(const Vector& x) const
{
	float rho = sqrt(x[0]*x[0] + x[2]*x[2]);
	if (rho - r >= abs(x[1]) - h/2.0) {
		return Vector(x[0]/rho, 0.0, x[2]/rho);
	} else {
		return Vector(0.0, sgn(x[1]), 0.0);
	}
}

float CylindricBoundary::volume() const
{
	return M_PI * r * r * h;
}

#ifdef HAVE_OPENGL
void CylindricBoundary::draw() const
{
	static GLUquadricObj* cylinder = gluNewQuadric();

	gluQuadricDrawStyle(cylinder, GLU_LINE);

	glPushMatrix();
	glRotated(90.0, 1.0, 0.0, 0.0);
	glTranslated(0.0, 0.0, -h/2.0);
	gluCylinder(cylinder, r, r, h, 50, 50);
	glPopMatrix();
}
#endif

float AnnulusBoundary::operator()(const Vector& x) const
{
	float rho = sqrt(x[0]*x[0] + x[2]*x[2]);
	return max<float>(max(rho - R, r - rho), abs(x[1]) - h/2.0);
}

Point AnnulusBoundary::get_random_position() const
{
	float x, y, z, rho;

	do {
		x = 2*R * (0.5 - drand48());
		z = 2*R * (0.5 - drand48());
		rho = sqrt(x*x + z*z);
	} while (rho <= r || rho >= R);

	y = h * (0.5 - drand48());

	return (1 - EPSILON) * Point(x, y, z);
}

Vector AnnulusBoundary::gradient(const Vector& x) const
{
	float rho = sqrt(x[0]*x[0] + x[2]*x[2]);
	if (max(rho - R, r - rho) >= abs(x[1]) - h/2.0) {
		Vector n(x[0]/rho, 0.0, x[2]/rho);
		return (rho - R > r - rho) ? n : -n;
	} else {
		return Vector(0.0, sgn(x[1]), 0.0);
	}
}

float AnnulusBoundary::volume() const
{
	return M_PI * (R*R - r*r) * h;
}

#ifdef HAVE_OPENGL
void AnnulusBoundary::draw() const
{
	static GLUquadricObj* cylinder = gluNewQuadric();

	gluQuadricDrawStyle(cylinder, GLU_LINE);

	glPushMatrix();
	glRotated(90.0, 1.0, 0.0, 0.0);
	glTranslated(0.0, 0.0, -h/2.0);
	gluCylinder(cylinder, r, r, h, 36, 36);
	gluCylinder(cylinder, R, R, h, 36, 36);
	glPopMatrix();
}
#endif

float AnnulusBoundary::predict_collision(const Particle& p) const
{
	float tolerance = 0.001;
	Vector x = p.position();
	Vector v = p.velocity();
	float g = p.growth_rate();
	float a, b, c, d, q, r0, R0, tr, tR, ty = t_stop;

	float tr1, tr2, tR1, tR2;
	float rho  = sqrt(x[0]*x[0] + x[2]*x[2]);
	float  vr  = (x[0]*v[0] + x[2]*v[2])/rho;

	/* compute escape time via top/bottom */

	if (v[1] + g >= 0.0)
		set_min<float>(ty, (h/2.0 - x[1] - g*t_curr - tolerance) / (v[1] + g));

	if (v[1] - g <= 0.0)
		set_min<float>(ty, (g*t_curr - h/2.0 - x[1] + tolerance) / (v[1] - g));

	if (ty < 0.0) ty = 0.0;

	/* compute escape time via cylinder */

	x[1] = 0.0; v[1] = 0.0; /* only use xz-plane */

	R0 = R - g*t_curr - tolerance;
	a = g*g - v*v;
	b = x*v + R0*g;
	c = R0*R0 - x*x;
	d = b*b - a*c;

	if (b > 0.0) {
		q = b + sqrt(d); tR1 = c/q; tR2 = q/a;
	} else {
		q = b - sqrt(d); tR1 = q/a; tR2 = c/q;
	}

	tR = tR1 > 0.0 ? tR1 : tR2;

	if (rho >= R0 || tR != tR)
		tR = vr + g > 0.0 ? 0.0 : INFINITY;

	if (abs(a) < EPSILON)
		tR = 0.5 * abs(c / b);

	r0 = r + g*t_curr + tolerance;
	b = r0*g - x*v;
	c = x*x - r0*r0;
	d = b*b + a*c;

	if (d > 0.0) {
		if (b > 0.0) {
			q = b + sqrt(d); tr1 = c/q; tr2 = -q/a;
		} else {
			q = b - sqrt(d); tr1 = -q/a; tr2 = c/q;
		}

		tr = tr1 > 0.0 ? tr1 : tr2;
	} else {
		tr = INFINITY; /* no intersection with inner cylinder */
	}

	if (rho <= r0)
		tr = vr - g < 0.0 ? 0.0 : INFINITY;

	if (abs(a) < EPSILON)
		tr = 0.5 * abs(c / b);

	if (tr < 0.0)
		tr = INFINITY; /* does not intersect inner cylinder */

	assert(ty  >= 0.0);
	assert(tr  >= 0.0);
	assert(tR  >= 0.0);
	assert(R - rho >= 0.5*tolerance);
	assert(rho - r >= 0.5*tolerance);

	return t_curr + min<float>(ty, tR, tr);
}

void AnnulusBoundary::collision_response(Particle& p) const
{
	float tolerance = 0.002;
	Vector x  = p.position();
	Vector v  = p.velocity();
	float  g  = p.growth_rate();
	float rho = sqrt(x[0]*x[0] + x[2]*x[2]);
	float  vr = (x[0]*v[0] + x[2]*v[2])/rho;

	if (h/2.0 - abs(x[1]) - g*t_curr < tolerance) {
		v[1] = abs(v[1]) > g ? -v[1] : -1.01 * sgn(x[1]) * g;
	}

	x[1] = 0.0; /* use only xz-plane position */

	if (R - rho - g*t_curr < tolerance) {
			v -= 2 * (vr >  g ? vr :  g) * x / rho;
	} else if (rho - r - g*t_curr < tolerance) {
			v -= 2 * (vr < -g ? vr : -g) * x / rho;
	}

	p.set_velocity(v);
}

float AnnulusWedgeBoundary::operator()(const Vector& x) const
{
	float rho = sqrt(x[0]*x[0] + x[2]*x[2]);
	float  d1 = max<float>(-x[2], -x*N); /* distance to wedge planes */
	float  d2 = max<float>(max(rho - R, r - rho), abs(x[1]) - h/2.0);
	return max<float>(d1, d2);
}

Point AnnulusWedgeBoundary::get_random_position() const
{
	float x, y, z, phi, rho;

	do {
		x = 2*R * (0.5 - drand48());
		z = R * drand48();
		rho = sqrt(x*x + z*z);
		phi = acos(x/rho);
	} while (rho <= r + 0.01 || rho >= R - 0.01
	        || phi < 0.035 || theta-phi < 0.035);

	y = 0.99 * h * (0.5 - drand48());

	return (1 - EPSILON) * Point(x, y, z);
}

Vector AnnulusWedgeBoundary::gradient(const Vector& x) const
{
	float rho = sqrt(x[0]*x[0] + x[2]*x[2]);
	float  d1 = max<float>(-x[2], -x*N); /* distance to wedge planes */
	float  d2 = max<float>(max(rho - R, r - rho), abs(x[1]) - h/2.0);

	if (d1 > d2) {
		if (x[2] < x*N)
			return Vector(0.0, 0.0, -1.0);
		else
			return N;
	} else if (max(rho - R, r - rho) >= abs(x[1]) - h/2.0) {
		Vector n(x[0]/rho, 0.0, x[2]/rho);
		return (rho - R > r - rho) ? n : -n;
	} else {
		return Vector(0.0, sgn(x[1]), 0.0);
	}
}

float AnnulusWedgeBoundary::volume() const
{
	return 0.5 * theta * (R*R - r*r) * h;
}

#ifdef HAVE_OPENGL
void AnnulusWedgeBoundary::draw() const
{
	static GLUquadricObj* cylinder = gluNewQuadric();

	gluQuadricDrawStyle(cylinder, GLU_LINE);

	glPushMatrix();
	glRotated(90.0, 1.0, 0.0, 0.0);
	glTranslated(0.0, 0.0, -h/2.0);
	gluCylinder(cylinder, r, r, h, 36, 36);
	gluCylinder(cylinder, R, R, h, 36, 36);
	glPopMatrix();

	glBegin(GL_LINES);
	for (int i = 0; i <= 10; i++) {
		glVertex3f(r, i/10.0*h-h/2.0, 0);
		glVertex3f(R, i/10.0*h-h/2.0, 0);

		glVertex3f(r*cos(theta), i/10.0*h-h/2.0, r*sin(theta));
		glVertex3f(R*cos(theta), i/10.0*h-h/2.0, R*sin(theta));

		glVertex3f(r + (i/10.0)*(R-r), -h/2.0, 0.0);
		glVertex3f(r + (i/10.0)*(R-r), +h/2.0, 0.0);

		float rr = r + (i/10.0)*(R-r);
		glVertex3f(rr*cos(theta), -h/2.0, rr*sin(theta));
		glVertex3f(rr*cos(theta), +h/2.0, rr*sin(theta));
	}
	glEnd();
}
#endif

float AnnulusWedgeBoundary::predict_collision(const Particle& p) const
{
	float tolerance = 0.001;
	Vector x = p.position();
	Vector v = p.velocity();
	float g = p.growth_rate();
	float a, b, c, d, q, r0, R0, tr, tR, tp = t_stop;

	float tr1, tr2, tR1, tR2;
	float rho  = sqrt(x[0]*x[0] + x[2]*x[2]);
	float  vr  = (x[0]*v[0] + x[2]*v[2])/rho;
	float  vw  = -v*N;

	/* compute escape time via top/bottom planes and wedge planes */

	if (v[1] + g >= -0.0) /* y = +h/2.0 plane */
		set_min<float>(tp, (h/2.0 - x[1] - g*t_curr - tolerance) / (v[1] + g));

	if (v[1] - g <= 0.0) /* y = -h/2.0 plane */
		set_min<float>(tp, (g*t_curr - h/2.0 - x[1] + tolerance) / (v[1] - g));

	if (v[2] - g <= 0.0) /* z = 0.0 plane */
		set_min<float>(tp, (x[2] - g*t_curr - tolerance) / (g - v[2]));

	if (vw + g >= 0.0) /* wedge plane at angle theta from z=0 plane */
		set_min<float>(tp, (x*N  - g*t_curr - tolerance) / (vw  +  g));

	if (tp < 0.0) tp = 0.0;

	/* compute escape time via cylinder */

	x[1] = 0.0; v[1] = 0.0; /* only use xz-plane */

	R0 = R - g*t_curr - tolerance;
	a = g*g - v*v;
	b = x*v + R0*g;
	c = R0*R0 - x*x;
	d = b*b - a*c;

	if (b > 0.0) {
		q = b + sqrt(d); tR1 = c/q; tR2 = q/a;
	} else {
		q = b - sqrt(d); tR1 = q/a; tR2 = c/q;
	}

	tR = tR1 > 0.0 ? tR1 : tR2;

	if (rho >= R0 || tR != tR)
		tR = vr + g > 0.0 ? 0.0 : INFINITY;

	if (abs(a) < EPSILON)
		tR = 0.5 * abs(c / b);

	r0 = r + g*t_curr + tolerance;
	b = r0*g - x*v;
	c = x*x - r0*r0;
	d = b*b + a*c;

	if (d > 0.0) {
		if (b > 0.0) {
			q = b + sqrt(d); tr1 = c/q; tr2 = -q/a;
		} else {
			q = b - sqrt(d); tr1 = -q/a; tr2 = c/q;
		}

		tr = tr1 > 0.0 ? tr1 : tr2;
	} else {
		tr = INFINITY; /* no intersection with inner cylinder */
	}

	if (rho <= r0)
		tr = vr - g < 0.0 ? 0.0 : INFINITY;

	if (abs(a) < EPSILON)
		tr = 0.5 * abs(c / b);

	if (tr < 0.0)
		tr = INFINITY; /* does not intersect inner cylinder */

	assert(tp  >= 0.0);
	assert(tr  >= 0.0);
	assert(tR  >= 0.0);
	assert(R - rho >= 0.5*tolerance);
	assert(rho - r >= 0.5*tolerance);

	return t_curr + min<float>(tp, tR, tr);
}

void AnnulusWedgeBoundary::collision_response(Particle& p) const
{
	float tolerance = 0.0012;
	Vector x  = p.position();
	Vector v  = p.velocity();
	float  g  = p.growth_rate();
	float rho = sqrt(x[0]*x[0] + x[2]*x[2]);
	float  vr = (x[0]*v[0] + x[2]*v[2])/rho;
	float  vw = -v*N;

	if (h/2.0 - abs(x[1]) - g*t_curr < tolerance) {
		v[1] = abs(v[1]) > g ? -v[1] : -1.01 * sgn(x[1]) * g;
	}

	if (x[2] - g*t_curr < tolerance && v[2] - g < 0.0) {
		v[2] = abs(v[2]) > g ? -v[2] : 1.01 * g;
	}

	if (x*N  - g*t_curr < tolerance && vw + g > 0.0) {
		v += 2 * (vw > g ? vw : g) * N;
	}

	x[1] = 0.0; /* use only xz-plane position */

	if (R - rho - g*t_curr < tolerance && vr > -g) {
			v -= 2 * (vr >  g ? vr :  g) * x / rho;
	} else if (rho - r - g*t_curr < tolerance && vr < g) {
			v -= 2 * (vr < -g ? vr : -g) * x / rho;
	}

	p.set_velocity(v);
}

float TorusBoundary::operator()(const Vector& x) const
{
	Vector p(x[0], 0.0, x[2]); p = R * p / norm(p);
	return norm(x - p) - r;
}

Vector TorusBoundary::gradient(const Vector& x) const
{
	Vector p(x[0], 0.0, x[2]); p = R * p / norm(p);
	return (x - p).normalized();
}

float TorusBoundary::predict_collision(const Particle& p) const
{
	Point  x = p.position();
	Vector v = p.velocity();
	Vector w = x - R * Point(x[0], 0.0, x[2]).normalized();

	float  g = p.growth_rate();
	float r0 = r - g * t_curr - 0.002;

	float a, b, c, d, q, t;

	a = v*v - g*g; b = v*w + g*r0; c = w*w - r0*r0; d = b*b - a*c;

	if (b > 0.0) {
		q = -(b + sqrt(d)); t = t_curr + c/q;
	} else {
		q = -(b - sqrt(d)); t = t_curr + q/a;
	}

	if (t != t) {
		t = t_curr; d = 1.0;
	}

	assert(norm(w) + g * t_curr < r); assert(d > -0.0);

	if (abs(a) < EPSILON) {
		return t_curr + 0.5 * abs(c / b);
	}

	return t;
}

float TorusBoundary::volume() const
{
	return 2 * M_PI * R * M_PI * r * r;
}

#ifdef HAVE_OPENGL
void TorusBoundary::draw() const
{
	glPushMatrix();
	glRotated(90.0, 1.0, 0.0, 0.0);
	glutWireTorus(r, R, 30, 60);
	glPopMatrix();
}
#endif
