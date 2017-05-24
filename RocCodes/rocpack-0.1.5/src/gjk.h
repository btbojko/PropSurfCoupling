#ifndef _GJK_H
#define _GJK_H

class Point;
class Vector;
class Particle;

bool contains(const Particle& A, const Point& x);
bool intersect(const Particle& A, const Particle& B, float t);
float distance(const Particle& A, const Particle& B, float t);
bool closest_points(const Particle& A, const Particle& B, Point& pA, Point& pB, float t);
void contact(const Particle& A, const Particle& B, Point& pA, Point& pB, Vector& N, float t);

#endif
