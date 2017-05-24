/* ROCSTAT - hotspot.c */

#include <stdio.h>

#include "amira.h"
#include "random.h"

enum { false, true };

extern FILE* output;

/* vertex connectivity */
int is_boundary(int x, int y, int z, int r){
	int dx, dy, dz;

	for(dz = -r; dz <= r; dz++){
		if(z + dz < 0 || z + dz >= mesh.size[2])
			return true;
		for(dy = -r; dy <= r; dy++){
			if(y + dy < 0 || y + dy >= mesh.size[1])
        return true;
			for(dx = -r; dx <= r; dx++){
				if(x + dx < 0 || x + dx >= mesh.size[0])
          return true;

				if((*mesh.voxel_content)(x+dx, y+dy, z+dz) == 0)
					return true;
			}
		}
	}
	return false;
}

void generate_hotspots(int n_hotspots, int r){
  int n;

  for (n = 0; n < n_hotspots; n++) {
			int x, y, z;
			do {
				x = (int) random_uniform(0, mesh.size[0]);
				y = (int) random_uniform(0, mesh.size[1]);
				z = (int) random_uniform(0, mesh.size[2]);
			} while(is_boundary(x,y,z,r));
			fprintf(output,"%4d %4d %4d %7.2f %7.2f %7.2f\n",
				x, y, z, x*mesh.psize[0], y*mesh.psize[1], z*mesh.psize[2]);
  }
}
