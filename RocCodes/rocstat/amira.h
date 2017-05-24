/* ROCSTAT - amira.h */

#ifndef _AMIRA_H
#define _AMIRA_H

typedef struct {
	float    bbox[3];
	float    psize[3];
	size_t   size[3];
	size_t   data_size;
	size_t   n_elements;
	size_t   n_particles;
	void*    data;
	float*   volume;
	unsigned char* bin_of;
	unsigned short (*voxel_content)(int,int,int);
} amira_mesh_t;

#define VOXEL_CONTENT_BYTE(x, y, z) \
	((unsigned char*) mesh.data)[ (z*mesh.size[1]+y)*mesh.size[0]+x ]

#define VOXEL_CONTENT_USHORT(x, y, z) \
	((unsigned short*) mesh.data)[ (z*mesh.size[1]+y)*mesh.size[0]+x ]

extern amira_mesh_t mesh;

unsigned short voxel_content_byte(int x, int y, int z);
unsigned short voxel_content_ushort(int x, int y, int z);

void load_amira_file(const char* filename, amira_mesh_t *m);

#endif
