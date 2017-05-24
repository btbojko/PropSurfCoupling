#ifndef _SCALAR_H
#define _SCALAR_H

#include "config.h"
#include "debug.h"

#include <cmath>
#include <cfloat>
#include <climits>
#include <cstdlib>

#ifndef EPSILON
#  define EPSILON FLT_EPSILON
#endif

#ifndef INFINITY
#  define INFINITY FLT_MAX
#endif

#define deg(x) ((x) * (180.0 / M_PI))
#define rad(x) ((x) / (180.0 / M_PI))

using std::abs;  /* use C++ version of abs()  */
using std::sqrt; /* use C++ version of sqrt() */

template <typename T>
inline T min(T a, T b)
{
	return a < b ? a : b;
}

template <typename T>
inline T max(T a, T b)
{
	return a > b ? a : b;
}

template <typename T>
inline T min(T a, T b, T c)
{
	return a < b ? min(a, c) : min(b, c);
}

template <typename T>
inline T max(T a, T b, T c)
{
	return a > b ? max(a, c) : max(b, c);
}

template <typename T>
inline void set_min(T& a, T b)
{
	if (a > b) a = b;
}

template <typename T>
inline void set_max(T& a, T b)
{
	if (a < b) a = b;
}

template <typename T>
inline T sgn(T x)
{
	return x > T(0) ? T(1) : x < T(0) ? T(-1) : T(0);
}

template <typename T>
inline bool inrange(T x, T min, T max)
{
	return min < max ? (x >= min) && (x <= max)
	                 : (x >= max) && (x <= min);
}

template <typename T>
void clamp(T& x, T min, T max)
{
	if (min < max) {
		if (x > max)
			x = max;
		else if (x < min)
			x = min;
	} else {
		if (x > min)
			x = min;
		else if (x < max)
			x = max;
	}
}

inline bool near_zero(float x, float epsilon = EPSILON)
{
	return abs(x) < epsilon;
}

inline bool nearly_equal(float x, float y, float epsilon = EPSILON)
{
	if (x == y) {
		return true;
	} else if (x * y == 0.0f) {
		float diff = abs(x - y);
		return diff < (epsilon * epsilon);
	} else {
		float diff = abs(x - y);
		float abs_x = abs(x), abs_y = abs(y);
		return diff / (abs_x + abs_y) < epsilon;
	}
}

#endif
