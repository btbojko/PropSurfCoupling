#ifndef _TRANSFORM_H
#define _TRANSFORM_H

#include "point.h"
#include "vector.h"
#include "matrix.h"
#include "quaternion.h"

class Transform {
  public:
	Transform() {}
	Transform(const float m[16]);
	Transform(const Transform& t)
		: basis(t.basis), origin(t.origin), type(t.type) {}
	Transform(const Matrix& M, const Point& x)
		: basis(M), origin(x), type(AFFINE) {}
	Transform(const Quaternion& q, const Point& x)
		: basis(q), origin(x), type(ROTATION | TRANSLATION) {}

	~Transform() {}

	Transform  operator*(const Transform& t);
	Transform& operator*=(const Transform& t);

	Point  operator*(const  Point& x) const;
	Point  operator()(const Point& x) const;

	/* vectors do not get translated */
	Vector operator*(const  Vector& x) const;
	Vector operator()(const Vector& x) const;

	/* coordinate transform for matrices */
	Matrix operator()(const Matrix& M) const;

	void set_identity();
	void translate(const Vector& v);
	void rotate(const Quaternion& q);
	void scale(const float s);
	void scale(const Vector& s);
	void scale(float x, float y, float z);

	const Point& get_origin() const;
	const Matrix& get_basis() const;

	void set_origin(const Point& p);
	void set_basis(const Matrix& M);
	void set_rotation(const Quaternion& q);

	void invert(const Transform& t);

  private:
	enum {
		IDENTITY    = 0x00,
		TRANSLATION = 0x01,
		ROTATION    = 0x02,
		SCALING     = 0x04,
		LINEAR      = ROTATION | SCALING,
		AFFINE      = LINEAR | TRANSLATION
	};

	Matrix basis; Point origin; unsigned int type;

	void mult(const Transform& t1, const Transform& t2);
	void mult_inverse_left(const Transform& t1, const Transform& t2);
	friend Transform inverse(const Transform& t);
};

Transform inverse(const Transform& t);

inline Transform& Transform::operator*=(const Transform& t)
{
	origin += basis * t.origin;
	basis *= t.basis;
	type |= t.type;
	return *this;
}

inline Transform Transform::operator*(const Transform& t)
{
	return Transform(basis * t.basis, (*this)(t.origin));
}

inline Point Transform::operator()(const Point& x) const
{
	return Point(basis*x + origin);
}

inline Point Transform::operator*(const Point& x) const
{
	return (*this)(x);
}

inline Vector Transform::operator()(const Vector& v) const
{
	return Vector(basis * v);
}

inline Vector Transform::operator*(const Vector& v) const
{
	return (*this)(v);
}

inline Matrix Transform::operator()(const Matrix& M) const
{
	Matrix MM = basis * M * (type & SCALING ? transpose(basis) : inverse(basis));
	return MM;
}

inline void Transform::set_identity()
{
	basis.set_identity(); origin.set_zero();
	type = IDENTITY;
}

inline void Transform::translate(const Vector& v)
{
	origin -= basis * v; type |= TRANSLATION;
}

inline void Transform::rotate(const Quaternion& q)
{
	if (norm(q) < EPSILON) return;
	basis *= Matrix(q); type |= ROTATION;
}

inline void Transform::scale(const float s)
{
	basis *= Matrix(s); type |= SCALING;
}

inline void Transform::scale(const Vector& v)
{
	basis *= Matrix(v); type |= SCALING;
}

inline void Transform::scale(float x, float y, float z)
{
	basis *= Matrix(x, y, z); type |= SCALING;
}

inline const Matrix& Transform::get_basis() const { return basis; }
inline const Point& Transform::get_origin() const { return origin; }

inline void Transform::set_basis(const Matrix& M) { basis = M; }
inline void Transform::set_origin(const Point& p) { origin = p; }
inline void Transform::set_rotation(const Quaternion& q) { basis.set_rotation(q); }

inline void Transform::invert(const Transform& t)
{
	basis = t.type & SCALING ? inverse(t.basis) : transpose(t.basis);
	origin = -(basis * t.origin); type = t.type;
}

inline void Transform::mult(const Transform& t1, const Transform& t2)
{
	basis = t1.basis * t2.basis; origin = t1(t2.origin);
	type = t1.type | t2.type;
}

inline void Transform::mult_inverse_left(const Transform& t1, const Transform& t2)
{
	Vector v = t2.origin - t1.origin;
	if (t1.type & SCALING) {
		Matrix inv = inverse(t1.basis); basis = inv * t2.basis; origin = inv * v;
	}
	else {
		basis = mult_transpose_left(t1.basis, t2.basis); origin = v * t1.basis;
	}
	type = t1.type | t2.type;
}

inline Transform inverse(const Transform& t)
{
	Matrix basis = t.type & Transform::SCALING ? inverse(t.basis) : transpose(t.basis);
	Point origin = -(basis * t.origin);
	return Transform(basis, origin);
}

#endif
