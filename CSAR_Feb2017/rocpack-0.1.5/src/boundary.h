#ifndef _BOUNDARY_H
#define _BOUNDARY_H

#include "vector.h"

class Particle;

class Boundary {
  public:
	virtual ~Boundary() {}

	virtual float operator()(const Vector& x) const = 0;

	/* provide default implementation for the grad */
	virtual Vector gradient(const Vector& x) const {
		const float h = 1e-6;
		static const Vector dx(h,0,0), dy(0,h,0), dz(0,0,h);

		float f  =  (*this)(x);
		float fx = ((*this)(x+dx) - f) / h;
		float fy = ((*this)(x+dy) - f) / h;
		float fz = ((*this)(x+dz) - f) / h;

		return Vector(fx, fy, fz);
	}

	virtual float predict_collision(const Particle& p) const;
	virtual void collision_response(Particle& p) const;
	virtual Point get_random_position() const;

	virtual bool is_periodic() const { return false; }
	virtual Vector dimensions() const = 0;

	virtual float volume() const = 0;

#ifdef HAVE_OPENGL
	virtual void draw() const = 0;
#endif
};

extern Vector period; /* for periodic boundaries */
extern Boundary* boundary;

/* types of boundary */

class PeriodicBoundary : public Boundary {
  public:
	PeriodicBoundary(float size) : size(Vector(size)) {}
	PeriodicBoundary(const Vector& size) : size(size) {}
	PeriodicBoundary(float width, float height, float depth)
		: size(Vector(width, height, depth)) {}

	float operator()(const Vector& x) const;
	Vector gradient(const Vector& x) const;

	float predict_collision(const Particle& p) const;
	void collision_response(Particle& p) const;

	virtual bool is_periodic() const { return true; }

	Vector dimensions() const { return size; }

	void reposition(Particle& p) const;

	float volume() const;
#ifdef HAVE_OPENGL
	void draw() const;
#endif
  private:
	Vector size;
};

class CubicBoundary : public Boundary {
  public:
	CubicBoundary(float a) : a(a) {}
	~CubicBoundary() {}

	float operator()(const Vector& x) const;
	Vector gradient(const Vector& x) const;

	float predict_collision(const Particle& p) const;
	void collision_response(Particle& p) const;

	Vector dimensions() const { return Vector(a, a, a); }
	float volume() const;
#ifdef HAVE_OPENGL
	void draw() const;
#endif

  private:
	float a;
};

class BoxBoundary : public Boundary {
  public:
	BoxBoundary(float size) : size(Vector(size)) {}
	BoxBoundary(const Vector& size) : size(size) {}
	BoxBoundary(float width, float height, float depth)
		: size(Vector(width, height, depth)) {}
	~BoxBoundary() {}

	float operator()(const Vector& x) const;
	Vector gradient(const Vector& x) const;

	float predict_collision(const Particle& p) const;
	void collision_response(Particle& p) const;
	Point get_random_position() const;

	Vector dimensions() const { return size; }

	float volume() const;
#ifdef HAVE_OPENGL
	void draw() const;
#endif

  private:
	Vector size;
};

class SphericBoundary : public Boundary {
  public:
	SphericBoundary(float r) : r(r) {}
	~SphericBoundary() {}

	float operator()(const Vector& x) const;
	Vector gradient(const Vector& x) const;

	Vector dimensions() const { return 2*Vector(r, r, r); }

	float volume() const;
#ifdef HAVE_OPENGL
	void draw() const;
#endif

  private:
	float r;
};

class CylindricBoundary : public Boundary {
  public:
	CylindricBoundary(float r, float h) : r(r), h(h) {}
	~CylindricBoundary() {}

	float operator()(const Vector& x) const;
	Vector gradient(const Vector& x) const;

	Vector dimensions() const { return Vector(2*r, h, 2*r); }

	float volume() const;
#ifdef HAVE_OPENGL
	void draw() const;
#endif

  private:
	float r, h;
};

class AnnulusBoundary : public Boundary {
  public:
	AnnulusBoundary(float r, float R, float h) : r(r), R(R), h(h) {}
	~AnnulusBoundary() {}

	float operator()(const Vector& x) const;
	Vector gradient(const Vector& x) const;

	void collision_response(Particle& p) const;
	float predict_collision(const Particle& p) const;
	Point get_random_position() const;

	Vector dimensions() const { return Vector(2*R, h, 2*R); }

	float volume() const;
#ifdef HAVE_OPENGL
	void draw() const;
#endif

  private:
	float r, R, h;
};

class AnnulusWedgeBoundary : public Boundary {
  public:
	AnnulusWedgeBoundary(float r, float R, float h, float theta)
		: r(r), R(R), h(h), theta(theta)
	{
		/* inward normal to wedge plane (other plane is z=0) */
		N = Vector(sin(theta), 0.0, -cos(theta));
	}
	~AnnulusWedgeBoundary() {}

	float operator()(const Vector& x) const;
	Vector gradient(const Vector& x) const;

	void collision_response(Particle& p) const;
	float predict_collision(const Particle& p) const;
	Point get_random_position() const;

	Vector dimensions() const { return Vector(2*R, h, 2*R); }

	float volume() const;
#ifdef HAVE_OPENGL
	void draw() const;
#endif

  private:
	float r, R, h, theta; Vector N;
};

class TorusBoundary : public Boundary {
  public:
	TorusBoundary(float R, float r) : R(R), r(r) {}
	~TorusBoundary() {}

	float operator()(const Vector& x) const;
	Vector gradient(const Vector& x) const;

	float predict_collision(const Particle& p) const;

	Vector dimensions() const { return Vector(2*(r+R), 2*r, 2*(R+r)); }

	float volume() const;
#ifdef HAVE_OPENGL
	void draw() const;
#endif

  private:
	float R, r;
};

#endif
