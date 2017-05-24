/* ROCSTAT - amira.c */

#include <math.h>
#include <stdio.h>
#include <errno.h>
#include <stdlib.h>
#include <libgen.h>
#include <string.h>

#include "amira.h"

amira_mesh_t mesh;
extern unsigned char particle_bins;

unsigned short (*voxel_content)(int,int,int);

unsigned short voxel_content_byte(int x, int y, int z){
  return (unsigned short)
		((unsigned char*)mesh.data)[(z*mesh.size[1]+y)*mesh.size[0]+x];
}

unsigned short voxel_content_ushort(int x, int y, int z){
  return ((unsigned short*)mesh.data)[(z*mesh.size[1]+y)*mesh.size[0]+x];
}

const char* find_and_jump(const char* buffer, const char* str){
	const char* found = strstr(buffer, str);
	return found ? found + strlen(str) : NULL;
}

void load_particle_data(const char* filename, amira_mesh_t *m){
	FILE* input; int i; char c, path[256], *dir, *file, *ext;

	strncpy(path, filename, sizeof(path) - 1);

	dir = dirname(path); file = basename(path);

	if(!(ext = strstr(file, ".am")) || strlen(ext) != 3){
		fprintf(stderr, "Input file must have the extension '.am'.\n");
		exit(1);
	}

	file[ strlen(file) - 3 ] = '\0';
	file = strcat(file, ".pdata");
	sprintf(path, "%s/%s", strdup(dir), strdup(file));

	if(!(input = fopen(path,"r")))
		return;

	m->n_particles = 0;
	while((c = getc(input)) != EOF)
		if(c == '\n') m->n_particles++;

	rewind(input);

	m->volume = (float*) malloc( m->n_particles * sizeof(float));
	m->bin_of = (unsigned char*) malloc( m->n_particles * sizeof(unsigned char));

	for(i = 0, m->volume[0] = 0; i < m->n_particles; i++){
		int index; float volume;

		if( fscanf(input, "%d %f", &index, &volume) != 2){
			fprintf(stderr, "%d\t\t%f\n", index, volume);
		}

		if(index < 0 || index > m->n_particles){
			fprintf(stderr, "Error reading %s: %d: index out of bounds.\n", path, index);
			exit(1);
		}

		m->volume[ index ] = volume;

		if(volume > m->volume[0])
			m->volume[0] = volume;
	}

	{ /* assign bins to each label */
		float r_max = pow(m->volume[0], 1.0/3.0);

		m->bin_of[0] = 0;

		for(i = 1; i < m->n_particles; i++){
			float r = pow(m->volume[i], 1.0/3.0);

			if(r != r_max)
				m->bin_of[i] = 1 + (unsigned char)(particle_bins * r/r_max);
			else
				m->bin_of[i] = particle_bins; /* biggest particle needs to be handled specially */
		}
	}

	fclose(input);
}

void load_amira_file(const char* filename, amira_mesh_t *m){
	FILE* input;
	char buffer[512];

	if(!(input = fopen(filename,"r"))){
		fprintf(stderr, "Error opening file %s: %s.\n", filename, strerror(errno));
		exit(errno);
	}

	if(fgets(buffer, 512, input) == NULL) {
		fprintf(stderr, "Error reading from file %s: %s.\n", filename, strerror(errno));
		exit(errno);
	}

	if (!strstr(buffer, "# AmiraMesh BINARY-LITTLE-ENDIAN 2.1")){
		fprintf(stderr, "Error: \'%s\' is not a proper Amira file.\n", filename);
		exit(EBADF);
	}

	while(fgets(buffer, 512, input) != NULL){
		const char *str = NULL;

		if((str = find_and_jump(buffer, "define Lattice"))){
			if(sscanf(str,"%lu %lu %lu", &m->size[0], &m->size[1], &m->size[2]) != 3){
				fprintf(stderr, "Error parsing input file.\n");
				exit(1);
			}
			m->n_elements = m->size[0] * m->size[1] * m->size[2];
		}

		if((str = find_and_jump(buffer, "Lattice {"))){
			char token[64]; strcpy(token, str);

			str = strtok(token, " \n");

			if(!strcmp(str, "byte")){
				m->data_size = sizeof(unsigned char);
				m->voxel_content = voxel_content_byte;
			} else if(!strcmp(str, "ushort")){
				m->data_size = sizeof(unsigned short);
				m->voxel_content = voxel_content_ushort;
			} else {
				fprintf(stderr, "Error: Unsupported data type: %s.\n", str);
				exit(1);
			}
		}

		if((str = find_and_jump(buffer, "BoundingBox"))){
			int i; float min[3], max[3];

			if(sscanf(str,"%f %f %f %f %f %f",
				&min[0], &max[0], &min[1], &max[1], &min[2], &max[2]) != 6){
				fprintf(stderr, "Error parsing input file.\n");
				exit(1);
			}

			for(i = 0; i < 3; i++){
				m->bbox[i] = (max[i]-min[i]);
				m->psize[i] = (m->bbox[i]/((float)(m->size[i])));
			}

#if 0
			/* only files with cubic voxels are allowed at the moment */
			if( fabs(m->psize[0] - m->psize[1]) > 1e-2 ||
			    fabs(m->psize[0] - m->psize[2]) > 1e-2){
				fprintf(stderr, "Error: Non-cubic voxels not yet implemented ");
				fprintf(stderr, "(%.4fum x %.4fum x %.4fum).\n",
					m->psize[0], m->psize[1], m->psize[2]);
				exit(1);
			}
#endif
		}

		/* break when reaching the end of the header section */
		if(strstr(buffer, "@1") == buffer)
			break;
	}

	m->data = malloc(m->n_elements * m->data_size);

	if(fread(m->data, m->data_size, m->n_elements, input) != m->n_elements || feof(input)){
		fprintf(stderr, "Error reading input file.\n");
		exit(errno);
	}

	fclose(input);

	if(m->data_size > 1)
		load_particle_data(filename, m);

}

