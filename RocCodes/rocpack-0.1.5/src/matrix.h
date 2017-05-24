#ifndef _MATRIX_H
#define _MATRIX_H

#include "vector.h"
#include "quaternion.h"

class Matrix {
  public:
	Matrix() { /* uninitialized */ }

	Matrix(float s) { set_scaling(s, s, s); }

	Matrix(const Vector& v){ set_scaling(v); }

	Matrix(float x, float y, float z) { set_scaling(x, y, z); }

	Matrix(const Quaternion& q) { set_rotation(q); }

	Matrix(const Vector& a, const Vector& b, const Vector& c)
	{
		set_value(a, b, c);
	}

	Matrix(__m128 rows[3])
	{
		data[0] = rows[0]; data[1] = rows[1]; data[2] = rows[2];
	}

	Matrix(__m128 row0, __m128 row1, __m128 row2)
	{
		data[0] = row0; data[1] = row1; data[2] = row2;
	}

	Matrix(float xx, float xy, float xz,
	       float yx, float yy, float yz,
	       float zx, float zy, float zz)
	{
		set_value(xx, xy, xz, yx, yy, yz, zx, zy, zz);
	}

	~Matrix() {}

	inline       Vector& operator[](int i)       { return (Vector&)(data[i]); }
	inline const Vector& operator[](int i) const { return (Vector&)(data[i]); }

	inline Matrix& mult_transpose_left(const Matrix& A)
	{
		_mm_multrl_ps(data, A.data); return *this;
	}

	inline Matrix& operator*=(const Matrix& A)
	{
		mult_transpose_left(A.transposed()); return *this;
	}

	inline void set_value(const Vector& a, const Vector& b, const Vector& c)
	{
		data[0] = a.get_value();
		data[1] = b.get_value();
		data[2] = c.get_value();
	}

	inline void set_value(float xx, float xy, float xz,
	                      float yx, float yy, float yz,
	                      float zx, float zy, float zz)
	{
		data[0] = _mm_set_ps(0, xz, xy, xx);
		data[1] = _mm_set_ps(0, yz, yy, yx);
		data[2] = _mm_set_ps(0, zz, zy, zx);
	}

	inline void set_identity()
	{
		set_value(1, 0, 0, 0, 1, 0, 0, 0, 1);
	}

	inline void set_scaling(float x, float y, float z)
	{
		set_value(x, 0, 0, 0, y, 0, 0, 0, z);
	}

	inline void set_scaling(const Vector& s)
	{
		set_value(s[0], 0, 0, 0, s[1], 0, 0, 0, s[2]);
	}

	inline void set_rotation(const Quaternion& q)
	{
		/* TODO: update to use SSE intrinsics */
		float d = q.norm2();
		float s = 2.0 / d;
		float x  = q[0], y  = q[1], z  = q[2], w = q[3];
		float xs = x*s,  ys = y*s,  zs = z*s;
		float x2 = x*xs, y2 = y*ys, z2 = z*zs;
		float xy = x*ys, xz = x*zs, yz = y*zs;
		float wx = w*xs, wy = w*ys, wz = w*zs;

		set_value( 1-(y2+z2),  (xy-wz) ,  (xz+wy),
					(xy+wz) , 1-(x2+z2),  (yz-wx),
					(xz-wy) ,  (yz+wx) , 1-(x2+y2));
	}

	inline float det() const
	{
		return _mm_cvtss_f32(_mm_dp_ps(data[0], _mm_cross3_ps(data[1], data[2]), 0x71));
	}

	inline Vector detv() const
	{
		return Vector(_mm_dp_ps(data[0], _mm_cross3_ps(data[1], data[2]), 0x7F));
	}

	inline void transpose() { _mm_transpose3_ps(); }

	inline Matrix transposed() const
	{
		Matrix M(data[0], data[1], data[2]);
		M.transpose();
		return M;
	}

	inline void invert() { _mm_invert3_ps(); }

	inline Matrix inverse() const
	{
		Matrix M(data[0], data[1], data[2]);
		M.invert();
		return M;
	}

  private:
	__m128 data[3]; /* rows */

	inline void _mm_multrl_ps(__m128 rows[3], const __m128 cols[3])
	{
		for (int i = 0; i < 3; i++) {
			__m128 x, y, z;
			x = _mm_dp_ps(rows[i], cols[0], 0x71);
			y = _mm_dp_ps(rows[i], cols[1], 0x71);
			z = _mm_dp_ps(rows[i], cols[2], 0x71);
			rows[i] = _mm_movelh_ps(_mm_unpacklo_ps(x, y), z);
		}
	}

	inline void _mm_transpose3_ps()
	{
		__m128 tmp0 = _mm_unpacklo_ps(data[0], data[1]);
		__m128 tmp2 = _mm_unpackhi_ps(data[0], data[1]);
		__m128 tmp3 = _mm_shuffle_ps(tmp2, data[2], _MM_SHUFFLE(2,3,3,2));
		__m128 tmp1 = _mm_unpacklo_ps(data[2], tmp3);
		__m128 tmp4 = _mm_unpackhi_ps(data[2], tmp3);
		data[0] = _mm_movelh_ps (tmp0, tmp1);
		data[1] = _mm_movehl_ps (tmp1, tmp0);
		data[2] = _mm_movelh_ps (tmp2, tmp4);
	}

	inline void _mm_invert3_ps()
	{
		__m128 tmp01 = _mm_shuffle_ps(data[0], data[0], _MM_SHUFFLE(3,0,2,1));
		__m128 tmp02 = _mm_shuffle_ps(data[0], data[0], _MM_SHUFFLE(3,1,0,2));
		__m128 tmp11 = _mm_shuffle_ps(data[1], data[1], _MM_SHUFFLE(3,0,2,1));
		__m128 tmp12 = _mm_shuffle_ps(data[1], data[1], _MM_SHUFFLE(3,1,0,2));
		__m128 tmp21 = _mm_shuffle_ps(data[2], data[2], _MM_SHUFFLE(3,0,2,1));
		__m128 tmp22 = _mm_shuffle_ps(data[2], data[2], _MM_SHUFFLE(3,1,0,2));

		__m128 col0 = _mm_sub_ps(_mm_mul_ps(tmp11,tmp22), _mm_mul_ps(tmp12,tmp21));
		__m128 col1 = _mm_sub_ps(_mm_mul_ps(tmp21,tmp02), _mm_mul_ps(tmp22,tmp01));
		__m128 col2 = _mm_sub_ps(_mm_mul_ps(tmp01,tmp12), _mm_mul_ps(tmp02,tmp11));

		__m128 det = _mm_dp_ps(data[0], col0, 0x7F);

		data[0] = _mm_div_ps(col0, det);
		data[1] = _mm_div_ps(col1, det);
		data[2] = _mm_div_ps(col2, det);

		_mm_transpose3_ps();
	}

	friend Matrix operator-(const Matrix& M);
	friend Matrix operator*(float s, const Matrix& A);
	friend Matrix operator*(const Matrix& M, float s);
	friend Matrix operator/(const Matrix& M, float s);
	friend Matrix operator+(const Matrix& A, const Matrix& B);
	friend Matrix operator-(const Matrix& A, const Matrix& B);
	friend Matrix operator*(const Matrix& A, const Matrix& B);
	friend Vector operator*(const Matrix& A, const Vector& v);
	friend Vector operator*(const Vector& v, const Matrix& A);

	friend Matrix mult_transpose_left(const Matrix& A, const Matrix& B);
};

inline float  det(const Matrix& A)  { return A.det(); }
inline Vector detv(const Matrix& A) { return A.detv(); }

inline Matrix transpose(const Matrix& A) { return A.transposed(); }
inline Matrix   inverse(const Matrix& A) { return A.inverse(); }

inline Matrix operator-(const Matrix& M)
{
	return Matrix(_mm_neg_ps(M.data[0]),
	              _mm_neg_ps(M.data[1]),
				  _mm_neg_ps(M.data[2]));
}

inline Matrix operator*(float s, const Matrix& M)
{
	__m128 sv = _mm_set1_ps(s);
	return Matrix(_mm_mul_ps(sv, M.data[0]),
	              _mm_mul_ps(sv, M.data[1]),
				  _mm_mul_ps(sv, M.data[2]));
}

inline Matrix operator*(const Matrix& M, float s)
{
	__m128 sv = _mm_set1_ps(s);
	return Matrix(_mm_mul_ps(sv, M.data[0]),
	              _mm_mul_ps(sv, M.data[1]),
				  _mm_mul_ps(sv, M.data[2]));
}

inline Matrix operator/(const Matrix& M, float s)
{
	__m128 sv = _mm_set1_ps(1.0 / s);
	return Matrix(_mm_mul_ps(sv, M.data[0]),
	              _mm_mul_ps(sv, M.data[1]),
				  _mm_mul_ps(sv, M.data[2]));
}

inline Matrix operator+(const Matrix& A, const Matrix& B)
{
	return Matrix(_mm_add_ps(A.data[0], B.data[0]),
	              _mm_add_ps(A.data[1], B.data[1]),
				  _mm_add_ps(A.data[2], B.data[2]));
}

inline Matrix operator-(const Matrix& A, const Matrix& B)
{
	return Matrix(_mm_sub_ps(A.data[0], B.data[0]),
	              _mm_sub_ps(A.data[1], B.data[1]),
				  _mm_sub_ps(A.data[2], B.data[2]));
}

inline Matrix mult_transpose_left(const Matrix& A, const Matrix& B)
{
	Matrix M(A.data[0], A.data[1], A.data[2]);
	return M.mult_transpose_left(B);
}

inline Matrix operator*(const Matrix& A, const Matrix& B)
{
	return mult_transpose_left(A, B.transposed());
}

inline Vector operator*(const Matrix& A, const Vector& v)
{
	__m128 x = _mm_dp_ps(A.data[0], v.get_value(), 0x71);
	__m128 y = _mm_dp_ps(A.data[1], v.get_value(), 0x71);
	__m128 z = _mm_dp_ps(A.data[2], v.get_value(), 0x71);
	return Vector(_mm_movelh_ps(_mm_unpacklo_ps(x, y), z));
}

inline Vector operator*(const Vector& v, const Matrix& A)
{
	return A.transposed() * v;
}

#endif
