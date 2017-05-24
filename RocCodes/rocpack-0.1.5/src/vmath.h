#ifndef _VMATH_H
#define _VMATH_H

#include <cmath>
#include <cfloat>

using std::abs;  /* Use C++ version of abs()  */
using std::sqrt; /* Use C++ version of sqrt() */

#if not defined(__SSE4_1__)
#  error "ERROR: Support for SSE4.1 SIMD instructions required, but not found"
#endif

#include <smmintrin.h> /* SSE4 intrinsics */

#define _mm_neg_ps(v) _mm_xor_ps(_mm_set1_ps(-0.0f), (v))
#define _mm_abs_ps(v) _mm_andnot_ps(_mm_set1_ps(-0.0f), (v))
#define _mm_norm3_ps(v) _mm_sqrt_ps(_mm_dp_ps((v), (v), 0x7F))
#define _mm_norm4_ps(v) _mm_sqrt_ps(_mm_dp_ps((v), (v), 0xFF))
#define _mm_madd_ps(v, a, b) _mm_add_ps(_mm_mul_ps((v), (a)), (b))

inline __m128 _mm_sgn_ps(__m128 v)
{
	__m128 cmp0 = _mm_cmplt_ps(v, _mm_setzero_ps());
	__m128 cmp1 = _mm_cmpgt_ps(v, _mm_setzero_ps());
	__m128 and0 = _mm_and_ps(cmp0, _mm_set1_ps(-1.0f));
	__m128 and1 = _mm_and_ps(cmp1, _mm_set1_ps(+1.0f));
	return _mm_or_ps(and0, and1);
}

/* vector cross product, for vectors as (x, y, z, 0) */

inline __m128 _mm_cross3_ps(__m128 a, __m128 b)
{
	__m128 tmp1 = _mm_shuffle_ps(a, a, _MM_SHUFFLE(3,0,2,1));
	__m128 tmp2 = _mm_shuffle_ps(b, b, _MM_SHUFFLE(3,1,0,2));
	__m128 tmp3 = _mm_shuffle_ps(a, a, _MM_SHUFFLE(3,1,0,2));
	__m128 tmp4 = _mm_shuffle_ps(b, b, _MM_SHUFFLE(3,0,2,1));
	return _mm_sub_ps(_mm_mul_ps(tmp1,tmp2), _mm_mul_ps(tmp3,tmp4));

	/* operations: 4 shuffles, 2 multiplications, 1 subtraction */
}

/* quaternion product of two quaternions (x, y, z, w) x (a, b, c, d) */

inline __m128 _mm_cross4_ps(__m128 xyzw, __m128 abcd)
{
	/* the product of the two quaternions is: */
	/* (X,Y,Z,W) = (xd+yc-zb+wa, -xc+yd+za+wb, xb-ya+zd+wc, -xa-yb-zc+wd) */

	__m128 wzyx = _mm_shuffle_ps(xyzw, xyzw, _MM_SHUFFLE(0,1,2,3));
	__m128 baba = _mm_shuffle_ps(abcd, abcd, _MM_SHUFFLE(0,1,0,1));
	__m128 dcdc = _mm_shuffle_ps(abcd, abcd, _MM_SHUFFLE(2,3,2,3));

	/* variable names below are for componens of result (X,Y,Z,W), nX for -X */

	/* znxwy  = (xb - ya, zb - wa, wd - zc, yd - xc) */
	__m128 ZnXWY = _mm_hsub_ps(_mm_mul_ps(xyzw, baba), _mm_mul_ps(wzyx, dcdc));

	/* xzynw  = (xd + yc, zd + wc, wb + za, yb + xa) */
	__m128 XZYnW = _mm_hadd_ps(_mm_mul_ps(xyzw, dcdc), _mm_mul_ps(wzyx, baba));

	/* _mm_shuffle_ps(XZYnW, ZnXWY, _MM_SHUFFLE(3,2,1,0)) = (xd + yc, zd + wc, wd - zc, yd - xc) */
	/* _mm_shuffle_ps(ZnXWY, XZYnW, _MM_SHUFFLE(2,3,0,1)) = (zb - wa, xb - ya, yb + xa, wb + za) */

	/* _mm_addsub_ps adds elements 1 and 3 and subtracts elements 0 and 2, so we get: */
	/* _mm_addsub_ps(*, *) = (xd+yc-zb+wa, xb-ya+zd+wc, wd-zc+yb+xa, yd-xc+wb+za)     */

	__m128 XZWY = _mm_addsub_ps(_mm_shuffle_ps(XZYnW, ZnXWY, _MM_SHUFFLE(3,2,1,0)),
	                            _mm_shuffle_ps(ZnXWY, XZYnW, _MM_SHUFFLE(2,3,0,1)));

	/* now we shuffle components in place and return the result */
	return _mm_shuffle_ps(XZWY, XZWY, _MM_SHUFFLE(2,1,3,0));

	/* operations: 6 shuffles, 4 multiplications, 3 compound additions/subtractions */
}

/* transpose a 4x4 matrix */

inline void _mm_transpose4_ps(__m128 row0, __m128 row1, __m128 row2, __m128 row3)
{
	__m128 tmp0 = _mm_unpacklo_ps(row0, row1);
	__m128 tmp1 = _mm_unpacklo_ps(row2, row3);
	__m128 tmp2 = _mm_unpackhi_ps(row0, row1);
	__m128 tmp3 = _mm_unpackhi_ps(row2, row3);
	row0 = _mm_movelh_ps (tmp0, tmp1); row1 = _mm_movehl_ps (tmp1, tmp0);
	row2 = _mm_movelh_ps (tmp2, tmp3); row3 = _mm_movehl_ps (tmp3, tmp2);

	/* operations: 4 unpacks, 4 moves */
}

/* 4x4 matrix times a 4 vector */

inline __m128 _mm_mat4mul_ps(__m128 row0, __m128 row1, __m128 row2, __m128 row3, __m128 v)
{
	__m128 x = _mm_dp_ps(row0, v, 0xF1); __m128 y = _mm_dp_ps(row1, v, 0xF1);
	__m128 z = _mm_dp_ps(row2, v, 0xF1); __m128 w = _mm_dp_ps(row3, v, 0xF1);
	return _mm_movelh_ps(_mm_unpacklo_ps(x,y), _mm_unpacklo_ps(z,w));
}

/* 3x3 matrix times a 3 vector + offset */

inline __m128 _mm_mat3madd_ps(__m128 row0, __m128 row1, __m128 row2, __m128 v, __m128 offset)
{
	__m128 x = _mm_dp_ps(row0, v, 0x71);
	__m128 y = _mm_dp_ps(row1, v, 0x71);
	__m128 z = _mm_dp_ps(row2, v, 0x71);
	return _mm_add_ps(_mm_movelh_ps(_mm_unpacklo_ps(x, y), z), offset);
}

#endif
