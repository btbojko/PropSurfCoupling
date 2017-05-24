#ifndef _COLLISION_H
#define _COLLISION_H

class Particle;

bool collision_check(const Particle& A, const Particle& B, float t);
float time_of_impact(const Particle& A, const Particle& B, float tmin, float tmax);
void collision_response(Particle& A, Particle& B);
void contact_network(Particle& p, FILE* out);

#endif
