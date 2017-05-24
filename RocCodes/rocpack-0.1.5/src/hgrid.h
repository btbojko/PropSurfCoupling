#ifndef _HGRID_H
#define _HGRID_H

#include <vector>
#include <unordered_map>

#include "scalar.h"

class Point;
class Particle;

typedef union {
	unsigned long int value;
	struct {
		short int x, y, z;
		unsigned short int level;
	} data;
} hash_t;

class HGrid {
  public:
	HGrid(float max_cell_size = 2.0, int max_level = 8);
	~HGrid() { grid.clear(); }

	void insert(Particle* p);
	void insert(Particle* p, hash_t hash);
	void erase(const Particle* p);
	void rehash(Particle* p);
	float next_rehash_time(const Particle* p) const;
	void find_neighbors(const Particle* p, std::vector<Particle*>& neighbors);

#ifdef HAVE_OPENGL
	void draw() const;
#endif

  private:
	const float MAX_CELL_SIZE;
	const unsigned int MAX_LEVEL;
	unsigned long int particles_at_level[32];
	unsigned int occupied_level_mask; /* bit mask of occupied levels */
	std::unordered_map<unsigned long int, Particle*> grid;

	hash_t make_hash(Particle* p);
	hash_t make_hash(const Point& p, float size);
};

extern HGrid* hgrid;

#endif
