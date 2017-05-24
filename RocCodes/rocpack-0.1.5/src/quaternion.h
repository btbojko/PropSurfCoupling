#ifndef _QUATERNION_H
#define _QUATERNION_H

#include "vector.h"

class Quaternion {
  public:
	inline Quaternion()
		: data(_mm_set_ps(1.0, 0.0, 0.0, 0.0)) {}

	inline Quaternion(float x, float y, float z, float w)
		: data(_mm_set_ps(w, z, y, x)) {}

	inline Quaternion(__m128 m) : data(m) {}

	inline Quaternion(const Vector& v)
		: data(_mm_set_ps(0.0, v[2], v[1], v[0])) {}

	inline Quaternion(const Vector& axis, float angle)
	{
		float s = sin(angle / 2.0) / axis.norm();
		data = _mm_set_ps(cos(angle / 2.0), axis[2]*s, axis[1]*s, axis[0]*s);
	}

	inline __m128 get_value() const { return data; }

	inline void set_value(float s) { data = _mm_set1_ps(s); }

	inline void set_value(float x, float y, float z, float w)
	{
		data = _mm_set_ps(w, z, y, x);
	}

	inline void set_identity() { data = _mm_set_ps(1.0, 0.0, 0.0, 0.0); }

	inline void set_rotation(const Vector& axis, float angle)
	{
		float s = sin(angle / 2.0) / axis.norm();
		data = _mm_set_ps(cos(angle / 2.0), axis[2]*s, axis[1]*s, axis[0]*s);
	}

	operator const float *() const { return (const float *) &data; }

	float  operator[](int k) const { return (&x)[k]; }
	float& operator[](int k)       { return (&x)[k]; }

	inline Quaternion& operator+=(float s)
	{
		data = _mm_add_ps(data, _mm_set1_ps(s)); return *this;
	}

	inline Quaternion& operator-=(float s)
	{
		data = _mm_sub_ps(data, _mm_set1_ps(s)); return *this;
	}

	inline Quaternion& operator*=(float s)
	{
		data = _mm_mul_ps(data, _mm_set1_ps(s)); return *this;
	}

	inline Quaternion& operator/=(float s)
	{
		data = _mm_div_ps(data, _mm_set1_ps(s)); return *this;
	}

	inline Quaternion& operator+=(const Quaternion& q)
	{
		data = _mm_add_ps(data, q.data); return *this;
	}

	inline Quaternion& operator-=(const Quaternion& q)
	{
		data = _mm_sub_ps(data, q.data); return *this;
	}

	inline Quaternion& operator^=(const Quaternion& q)
	{
		data = _mm_cross4_ps(data, q.data); return *this;
	}

	inline bool near_zero() const { return norm() < EPSILON; }

	inline float norm() const
	{
		return _mm_cvtss_f32(_mm_sqrt_ss(_mm_dp_ps(data, data, 0xF1)));
	}

	inline float rnorm() const
	{
		return _mm_cvtss_f32(_mm_rsqrt_ss(_mm_dp_ps(data, data, 0xF1)));
	}

	inline float norm2() const
	{
		return _mm_cvtss_f32(_mm_dp_ps(data, data, 0xF1));
	}

	inline Quaternion normv() const
	{
		return Quaternion(_mm_sqrt_ps(_mm_dp_ps(data, data, 0xFF)));
	}

	inline Quaternion rnormv() const
	{
		return Quaternion(_mm_rsqrt_ps(_mm_dp_ps(data, data, 0xFF)));
	}

	inline Quaternion norm2v() const
	{
		return Quaternion(_mm_dp_ps(data, data, 0xFF));
	}

	inline void normalize()
	{
		data = _mm_div_ps(data, _mm_norm4_ps(data));
	}

	inline Quaternion normalized() const
	{
		return Quaternion(_mm_div_ps(data, _mm_norm4_ps(data)));
	}

	inline void conjugate()
	{
		static const __attribute__ ((aligned(16)))
			unsigned int _mask[4] = {0x80000000,0x80000000,0x80000000,0};
		data = _mm_xor_ps(data, _mm_load_ps((const float *)_mask));
	}

	inline Quaternion conjugated() const
	{
		static const __attribute__ ((aligned(16)))
			unsigned int _mask[4] = {0x80000000,0x80000000,0x80000000,0};
		return Quaternion(_mm_xor_ps(data, _mm_load_ps((const float *)_mask)));
	}

	inline void invert()
	{
		conjugate(); data = _mm_div_ps(data, _mm_dp_ps(data, data, 0xFF));
	}

	inline Quaternion inverse() const
	{
		Quaternion q(data); q.invert(); return q;
	}

	inline int min_element() const
	{
		return x < w ? (x < z ? (x < y ? 0 : 1) : (y < z ? 1 : 2)) :
		               (w < z ? (w < y ? 3 : 1) : (y < z ? 1 : 2));
	}

	inline int max_element() const
	{
		return x > w ? (x > z ? (x > y ? 0 : 1) : (y > z ? 1 : 2)) :
		               (w > z ? (w > y ? 0 : 1) : (y > z ? 1 : 2));
	}

	static Quaternion random();
	static Quaternion identity() { return Quaternion(0.0, 0.0, 0.0, 1.0); }

  protected:
	union {
		__m128 data;
		struct { float x, y, z, w; };
	} __attribute__ ((aligned(16)));
};

inline float norm(const Quaternion& q) { return q.norm(); }
inline float norm2(const Quaternion& q) { return q.norm2(); }
inline float rnorm(const Quaternion& q) { return q.rnorm(); }

inline Quaternion  normv(const Quaternion& q) { return  q.normv(); }
inline Quaternion rnormv(const Quaternion& q) { return q.rnormv(); }
inline Quaternion norm2v(const Quaternion& q) { return q.norm2v(); }

inline Quaternion inverse(const Quaternion& q) { return q.inverse(); }
inline Quaternion conjugate(const Quaternion& q) { return q.conjugated(); }

inline bool near_zero(const Quaternion& q) { return q.near_zero(); }

inline Quaternion operator-(const Quaternion& q)
{
	return Quaternion(_mm_neg_ps(q.get_value()));
}

inline Quaternion operator+(float s, const Quaternion& q)
{
	return Quaternion(_mm_add_ps(q.get_value(), _mm_set1_ps(s)));
}

inline Quaternion operator+(const Quaternion& q, float s)
{
	return Quaternion(_mm_add_ps(q.get_value(), _mm_set1_ps(s)));
}

inline Quaternion operator-(float s, const Quaternion& q)
{
	return Quaternion(_mm_sub_ps(_mm_set1_ps(s), q.get_value()));
}

inline Quaternion operator-(const Quaternion& q, float s)
{
	return Quaternion(_mm_sub_ps(q.get_value(), _mm_set1_ps(s)));
}

inline Quaternion operator*(float s, const Quaternion& q)
{
	return Quaternion(_mm_mul_ps(q.get_value(), _mm_set1_ps(s)));
}

inline Quaternion operator*(const Quaternion& q, float s)
{
	return Quaternion(_mm_mul_ps(q.get_value(), _mm_set1_ps(s)));
}

inline Quaternion operator/(const Quaternion& q, float s)
{
	return Quaternion(_mm_div_ps(q.get_value(), _mm_set1_ps(s)));
}

inline Quaternion operator+(const Quaternion& q1, const Quaternion& q2)
{
	return Quaternion(_mm_add_ps(q1.get_value(), q2.get_value()));
}

inline Quaternion operator-(const Quaternion& q1, const Quaternion& q2)
{
	return Quaternion(_mm_sub_ps(q1.get_value(), q2.get_value()));
}

inline float dot(const Quaternion& q1, const Quaternion& q2)
{
	return _mm_cvtss_f32(_mm_dp_ps(q1.get_value(), q2.get_value(), 0xF1));
}

inline Quaternion dotv(const Quaternion& q1, const Quaternion& q2)
{
	return Quaternion(_mm_dp_ps(q1.get_value(), q2.get_value(), 0xFF));
}

inline Quaternion cross(const Quaternion& q1, const Quaternion& q2)
{
	return Quaternion(_mm_cross4_ps(q1.get_value(), q2.get_value()));
}

inline float operator*(const Quaternion& q1, const Quaternion& q2)
{
	return _mm_cvtss_f32(_mm_dp_ps(q1.get_value(), q2.get_value(), 0xF1));
}

inline Quaternion operator^(const Quaternion& q1, const Quaternion& q2)
{
	return Quaternion(_mm_cross4_ps(q1.get_value(), q2.get_value()));
}

inline Vector operator^(const Vector& v, const Quaternion& q)
{
	return Vector(_mm_cross4_ps(v.get_value(), q.get_value()));
}

inline Vector operator^(const Quaternion& q, const Vector& v)
{
	return Vector(_mm_cross4_ps(q.get_value(), v.get_value()));
}

/* multiply and divide element by element */

inline Quaternion elem_mul(const Quaternion& q1, const Quaternion& q2)
{
	return Quaternion(_mm_mul_ps(q1.get_value(), q2.get_value()));
}

inline Quaternion elem_div(const Quaternion& q1, const Quaternion& q2)
{
	return Quaternion(_mm_div_ps(q1.get_value(), q2.get_value()));
}

inline Quaternion Quaternion::random() {
    double x0 = drand48();
    double r1 = sqrt(1.0 - x0), r2 = sqrt(x0);
    double t1 = 2.0 * M_PI * drand48(), t2 = 2.0 * M_PI * drand48();
    double c1 = cos(t1), s1 = sin(t1);
    double c2 = cos(t2), s2 = sin(t2);
    return Quaternion(s1 * r1, c1 * r1, s2 * r2, c2 * r2);
}

#endif /* _QUATERNION_H_ */
