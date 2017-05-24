#ifndef _VECTOR_H
#define _VECTOR_H

#include "vmath.h"
#include "scalar.h"

class Point;

class Vector {
  public:
	inline Vector()
		: data(_mm_setzero_ps()) {}

	inline Vector(__m128 m)
		: data(m) {}

	inline Vector(float s)
		: data(_mm_set1_ps(s)) {}

	inline Vector(float x, float y, float z)
		: data(_mm_set_ps(0, z, y, x)) {}

	__m128 get_value() const { return data; }

	inline void set_zero() { data = _mm_setzero_ps(); }

	inline void set_value(float s)
	{
		data = _mm_set1_ps(s);
	}

	inline void set_value(__m128 xyzw)
	{
		data = xyzw;
	}

	inline void set_value(float x, float y, float z)
	{
		data = _mm_set_ps(0, z, y, x);
	}

	operator const float *() const { return (const float *) &data; }

	float  operator[](int k) const { return (&x)[k]; }
	float& operator[](int k)       { return (&x)[k]; }

	inline Vector& operator+=(float s)
	{
		data = _mm_add_ps(data, _mm_set1_ps(s)); return *this;
	}

	inline Vector& operator-=(float s)
	{
		data = _mm_sub_ps(data, _mm_set1_ps(s)); return *this;
	}

	inline Vector& operator*=(float s)
	{
		data = _mm_mul_ps(data, _mm_set1_ps(s)); return *this;
	}

	inline Vector& operator/=(float s)
	{
		data = _mm_div_ps(data, _mm_set1_ps(s)); return *this;
	}


	inline Vector& operator+=(const Vector& v)
	{
		data = _mm_add_ps(data, v.data); return *this;
	}

	inline Vector& operator-=(const Vector& v)
	{
		data = _mm_sub_ps(data, v.data); return *this;
	}


	inline Vector& elem_mul(const Vector& v)
	{
		data = _mm_mul_ps(data, v.data); return *this;
	}

	inline Vector& elem_div(const Vector& v)
	{
		data = _mm_div_ps(data, v.data); return *this;
	}


	inline float norm() const
	{
		return _mm_cvtss_f32(_mm_sqrt_ss(_mm_dp_ps(data, data, 0x71)));
	}

	inline float rnorm() const
	{
		return _mm_cvtss_f32(_mm_rsqrt_ss(_mm_dp_ps(data, data, 0x71)));
	}

	inline float norm2() const
	{
		return _mm_cvtss_f32(_mm_dp_ps(data, data, 0x71));
	}


	inline Vector normv() const
	{
		return Vector(_mm_norm3_ps(data));
	}

	inline Vector rnormv() const
	{
		return Vector(_mm_rcp_ps(_mm_norm3_ps(data)));
	}

	inline Vector norm2v() const
	{
		return Vector(_mm_dp_ps(data, data, 0x7F));
	}


	inline void normalize()
	{
		data = _mm_div_ps(data, _mm_norm3_ps(data));
	}

	inline Vector normalized() const
	{
		return Vector(_mm_div_ps(data, _mm_norm3_ps(data)));
	}


	inline bool near_zero() const
	{
		return norm() < EPSILON;
	}


	inline Vector abs() const
	{
		return Vector(_mm_abs_ps(data));
	}

	inline Vector sgn() const
	{
		return Vector(_mm_sgn_ps(data));
	}


	inline float min() const
	{
		return x < z ? (x < y ? x : y) : (y < z ? y : z);
	}

	inline float max() const
	{
		return x > z ? (x > y ? x : y) : (y > z ? y : z);
	}


	inline int min_element() const
	{
		return x < z ? (x < y ? 0 : 1) : (y < z ? 1 : 2);
	}

	inline int max_element() const
	{
		return x > z ? (x > y ? 0 : 1) : (y > z ? 1 : 2);
	}

	inline int closest_axis() const
	{
		return Vector(_mm_abs_ps(data)).max_element();
	}

	static Vector random() /* unit vector in random direction */
	{
		float z = 2.0 * drand48() - 1.0;
		float r = sqrt(1.0 - z * z);
		float t = 2.0 * M_PI * drand48();
		return Vector(r * cos(t), r * sin(t), z);
	}

  protected:
	union {
		__m128 data;
		struct { float x, y, z, w; };
	} __attribute__ ((aligned(16)));
};

inline float norm(const Vector& v) { return v.norm(); }
inline float norm2(const Vector& v) { return v.norm2(); }
inline float rnorm(const Vector& v) { return v.rnorm(); }

inline Vector normv(const Vector& v) { return v.normv(); }
inline Vector rnormv(const Vector& v) { return v.rnormv(); }
inline Vector norm2v(const Vector& v) { return v.norm2v(); }

inline bool near_zero(const Vector& v) { return v.near_zero(); }

/* other common operations */
inline float  min(const Vector& v) { return v.min(); }
inline float  max(const Vector& v) { return v.max(); }
inline Vector abs(const Vector& v) { return v.abs(); }
inline Vector sgn(const Vector& v) { return v.sgn(); }

inline Vector operator-(const Vector& v)
{
	return Vector(_mm_neg_ps(v.get_value()));
}

inline Vector operator+(float s, const Vector& v)
{
	return Vector(_mm_add_ps(v.get_value(), _mm_set1_ps(s)));
}

inline Vector operator+(const Vector& v, float s)
{
	return Vector(_mm_add_ps(v.get_value(), _mm_set1_ps(s)));
}

inline Vector operator-(float s, const Vector& v)
{
	return Vector(_mm_sub_ps(v.get_value(), _mm_set1_ps(s)));
}

inline Vector operator-(const Vector& v, float s)
{
	return Vector(_mm_sub_ps(v.get_value(), _mm_set1_ps(s)));
}

inline Vector operator*(float s, const Vector& v)
{
	return Vector(_mm_mul_ps(v.get_value(), _mm_set1_ps(s)));
}

inline Vector operator*(const Vector& v, float s)
{
	return Vector(_mm_mul_ps(v.get_value(), _mm_set1_ps(s)));
}

inline Vector operator/(const Vector& v, float s)
{
	return Vector(_mm_div_ps(v.get_value(), _mm_set1_ps(s)));
}

inline Vector operator+(const Vector& v1, const Vector& v2)
{
	return Vector(_mm_add_ps(v1.get_value(), v2.get_value()));
}

inline Vector operator-(const Vector& v1, const Vector& v2)
{
	return Vector(_mm_sub_ps(v1.get_value(), v2.get_value()));
}

inline float  operator*(const Vector& v1, const Vector& v2)
{
	return _mm_cvtss_f32(_mm_dp_ps(v1.get_value(), v2.get_value(), 0x71));
}

inline Vector operator^(const Vector& v1, const Vector& v2)
{
	return Vector(_mm_cross3_ps(v1.get_value(), v2.get_value()));
}

inline float dot(const Vector& v1, const Vector& v2)
{
	return _mm_cvtss_f32(_mm_dp_ps(v1.get_value(), v2.get_value(), 0x71));
}

inline Vector dotv(const Vector& v1, const Vector& v2)
{
	return Vector(_mm_dp_ps(v1.get_value(), v2.get_value(), 0x7F));
}

inline Vector cross(const Vector& v1, const Vector& v2)
{
	return Vector(_mm_cross3_ps(v1.get_value(), v2.get_value()));
}

/* multiply and divide element by element */

inline Vector elem_mul(const Vector& v1, const Vector& v2)
{
	return Vector(_mm_mul_ps(v1.get_value(), v2.get_value()));
}

inline Vector elem_div(const Vector& v1, const Vector& v2)
{
	return Vector(_mm_div_ps(v1.get_value(), v2.get_value()));
}

inline Vector min(const Vector& v1, const Vector& v2)
{
	return Vector(_mm_min_ps(v1.get_value(), v2.get_value()));
}

inline Vector max(const Vector& v1, const Vector& v2)
{
	return Vector(_mm_max_ps(v1.get_value(), v2.get_value()));
}

inline float angle(const Vector& v1, const Vector& v2)
{
    return acos((v1*v2) / sqrt(norm2(v1) * norm2(v2)));
}

#endif
