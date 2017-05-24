/* ROCSTAT - stat.c */

#include <omp.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* for memset */

#include "amira.h"
#include "random.h"

extern FILE*         output;
extern double        angle;
extern double        max_distance;
extern unsigned int  bins;
extern unsigned int  points_per_bin;
extern unsigned char particle_bins;
extern char*         slice;

void calculate_Sij(){
	register int i, j, k;
	unsigned int ***Sij;
	float step =  max_distance / ((bins-1) * mesh.psize[0]);

	particle_bins = 1; /* two phase pack by default */

	switch (mesh.data_size) {
		case 1:
			for(i = 0; i < mesh.n_elements; i++) {
				int tag = ((unsigned char*)mesh.data)[i];
				if (tag > particle_bins)
					particle_bins = tag;
			}
			break;
		case 2:
			for(i = 0; i < mesh.n_elements; i++) {
				int tag = ((unsigned short*)mesh.data)[i];
				if (tag > particle_bins)
					particle_bins = tag;
			}
			break;
	}

	Sij = (unsigned int***) malloc((particle_bins+1) * sizeof(unsigned int**));

	for(i = 0; i <= particle_bins; i++){
		Sij[i] = (unsigned int**) malloc((particle_bins+1) * sizeof(unsigned int*));

		for(j = 0; j <= particle_bins; j++){
			Sij[i][j] = (unsigned int*) malloc(bins * sizeof(unsigned int));
			memset((void*)Sij[i][j], 0, bins * sizeof(unsigned int));
		}
	}

#pragma omp parallel for default(shared) private(i,j)
	for(i = 0; i < bins; i++){
		unsigned short st[3];
		double distance = i * step;
		for(j = 0; j < points_per_bin; j++){
			int bin_x, bin_xdx; double x[3], v[3];

			/* generate a point inside our volume */
			x[0] = random_uniform_r(st, 0, mesh.size[0]);
			x[1] = random_uniform_r(st, 0, mesh.size[1]);
			x[2] = random_uniform_r(st, 0, mesh.size[2]);

			/* save contents of this location */
			bin_x = (*mesh.voxel_content)((int)x[0], (int)x[1], (int)x[2]);

			random_unit_vector_r(st, v);

			/* add random vector to our position */
			x[0] += distance * v[0];
			x[1] += distance * v[1];
			x[2] += distance * v[2];

			/* if x+dx is outside our volume, discard it */
			if(x[0] < 0 || x[0] >= mesh.size[0] ||
			   x[1] < 0 || x[1] >= mesh.size[1] ||
			   x[2] < 0 || x[2] >= mesh.size[2]){
				j = j-1;
				continue;
			}

			/* save contents of second location */
			bin_xdx = (*mesh.voxel_content)((int)x[0], (int)x[1], (int)x[2]);

#pragma omp atomic
			/* increase the right bin */
			Sij[bin_x][bin_xdx][i]++;
		}
	}

	/* write output */
	for(k = 0; k < bins; k++){
		fprintf(output, "%8.3f ", max_distance * k / (bins-1));
		for(i = 0; i <= particle_bins; i++)
			for(j = 0; j <= particle_bins; j++)
				fprintf(output,"%.6f ", (double)(Sij[i][j][k])/points_per_bin);
		fprintf(output,"\n");
	}

	for(i = 0; i < particle_bins; i++){
		for(j = 0; j < particle_bins; j++)
			free(Sij[i][j]);

		free(Sij[i]);
	}
}

void calculate_Sijk(){
	register int i, j, k;
	unsigned int **Sijk[8];
	double step =  max_distance / ((bins-1) * mesh.psize[0]);

	for(i = 0; i < 8; i++){
		Sijk[i] = (unsigned int**) malloc(bins * sizeof(unsigned int**));

		for(j = 0; j < bins; j++){
			Sijk[i][j] = (unsigned int*) malloc(bins * sizeof(unsigned int));
			memset((void*)Sijk[i][j], 0, bins * sizeof(unsigned int));
		}
	}

#pragma omp parallel for default(shared) private(i,j,k)
	for(i = 0; i < bins; i++){
		unsigned short st[3];
		double r1 = i * step;
		for(j = 0; j < bins; j++){
			double r2 = j * step;
			for(k = 0; k < points_per_bin; k++){
				int bin, content_i, content_j, content_k;
				double x[3], v1[3], v2[3];

				/* generate a point inside our volume */
				x[0] = random_uniform_r(st, 0, mesh.size[0]);
				x[1] = random_uniform_r(st, 0, mesh.size[1]);
				x[2] = random_uniform_r(st, 0, mesh.size[2]);

				/* save contents of this location */
				content_i = (*mesh.voxel_content)((int)x[0], (int)x[1], (int)x[2]) > 0;

				/* generate two random orthogonal vectors */
				random_unit_vector_r(st, v1);
				random_normal_r(st, v2, v1);

				/* change v2 as to be at an angle alpha with v1 */
				v2[0] = cos(angle) * v1[0] + sin(angle) * v2[0];
				v2[1] = cos(angle) * v1[1] + sin(angle) * v2[1];
				v2[2] = cos(angle) * v1[2] + sin(angle) * v2[2];

				/* add vectors to our position */
				v1[0] = x[0] + r1 * v1[0];
				v1[1] = x[1] + r1 * v1[1];
				v1[2] = x[2] + r1 * v1[2];

				/* if x+v1 is outside our volume, discard it */
				if(v1[0] < 0 || v1[0] >= mesh.size[0] ||
				   v1[1] < 0 || v1[1] >= mesh.size[1] ||
				   v1[2] < 0 || v1[2] >= mesh.size[2] ){
					k = k-1;
					continue;
				}

				v2[0] = x[0] + r2 * v2[0];
				v2[1] = x[1] + r2 * v2[1];
				v2[2] = x[2] + r2 * v2[2];

				/* if x+v2 is outside our volume, discard it */
				if(v2[0] < 0 || v2[0] >= mesh.size[0] ||
				   v2[1] < 0 || v2[1] >= mesh.size[1] ||
				   v2[2] < 0 || v2[2] >= mesh.size[2] ){
					k = k-1;
					continue;
				}

				/* save contents of other locations */
				content_j = (*mesh.voxel_content)((int)v1[0], (int)v1[1], (int)v1[2]) > 0;
				content_k = (*mesh.voxel_content)((int)v2[0], (int)v2[1], (int)v2[2]) > 0;

				bin = 2*content_i + 4*content_j + content_k;

#pragma omp atomic
				/* increase the right bin */
				Sijk[bin][i][j]++;
			}
		}
	}

	/* write result to output file */
	for(i = 0; i < bins; i++){
		for(j = 0; j < bins; j++){
			fprintf(output, "%7.3f %7.3f", i*step*mesh.psize[0], j*step*mesh.psize[0]);
			for(k = 0; k < 8; k++)
				fprintf(output," %f", (double) (Sijk[k][i][j])/points_per_bin);
			fprintf(output, "\n");
		}
		fprintf(output, "\n");
	}

	for(i = 0; i < 8; i++){
		for(j = 0; j < bins; j++)
			free(Sijk[i][j]);

		free(Sijk[i]);
	}
}

void calculate_Sij_bins(){
	register int i, j, k;
	unsigned int ***Sij;
	float step =  max_distance / ((bins-1) * mesh.psize[0]);

	if(mesh.bin_of == NULL){
		fprintf(stderr, "Particle data not found.\n");
		exit(1);
	}

	Sij = (unsigned int***) malloc((particle_bins+1) * sizeof(unsigned int**));

	for(i = 0; i <= particle_bins; i++){
		Sij[i] = (unsigned int**) malloc((particle_bins+1) * sizeof(unsigned int*));

		for(j = 0; j <= particle_bins; j++){
			Sij[i][j] = (unsigned int*) malloc(bins * sizeof(unsigned int));
			memset((void*)Sij[i][j], 0, bins * sizeof(unsigned int));
		}
	}

#pragma omp parallel for default(shared) private(i,j)
	for(i = 0; i < bins; i++){
		unsigned short st[3];
		double distance = i * step;
		for(j = 0; j < points_per_bin; j++){
			int content_x, bin_x, content_xdx, bin_xdx; double x[3], v[3];

			/* generate a point inside our volume */
			x[0] = random_uniform_r(st, 0, mesh.size[0]);
			x[1] = random_uniform_r(st, 0, mesh.size[1]);
			x[2] = random_uniform_r(st, 0, mesh.size[2]);

			/* save contents of this location */
			content_x = (*mesh.voxel_content)((int)x[0], (int)x[1], (int)x[2]);

			bin_x = mesh.bin_of[ content_x ];

			random_unit_vector_r(st, v);

			/* add random vector to our position */
			x[0] += distance * v[0];
			x[1] += distance * v[1];
			x[2] += distance * v[2];

			/* if x+dx is outside our volume, discard it */
			if(x[0] < 0 || x[0] >= mesh.size[0] ||
			   x[1] < 0 || x[1] >= mesh.size[1] ||
			   x[2] < 0 || x[2] >= mesh.size[2]){
				j = j-1;
				continue;
			}

			/* save contents of second location */
			content_xdx = (*mesh.voxel_content)((int)x[0], (int)x[1], (int)x[2]);

			bin_xdx = mesh.bin_of[ content_xdx ];

#pragma omp atomic
			/* increase the right bin */
			Sij[bin_x][bin_xdx][i]++;
		}
	}

	/* write output */
	for(k = 0; k < bins; k++){
		fprintf(output, "%8.3f ", max_distance * k / (bins-1));
		for(i = 0; i <= particle_bins; i++)
			for(j = 0; j <= particle_bins; j++)
				fprintf(output,"%.6f ", (double)(Sij[i][j][k])/points_per_bin);
		fprintf(output,"\n");
	}

	for(i = 0; i < particle_bins; i++){
		for(j = 0; j < particle_bins; j++)
			free(Sij[i][j]);

		free(Sij[i]);
	}
}

void calculate_Sijk_angle(){
	register int i, j, k;
	unsigned int **Sijk[8];
	double step =  max_distance / (bins * mesh.psize[0]);

	for(i = 0; i < 8; i++){
		Sijk[i] = (unsigned int**) malloc(bins * sizeof(unsigned int**));

		for(j = 0; j <= bins; j++) {
			Sijk[i][j] = (unsigned int*) malloc(181 * sizeof(unsigned int));
			memset((void*)Sijk[i][j], 0, 181 * sizeof(unsigned int));
		}
	}

#pragma omp parallel for default(shared) private(i,j,k)
	for(i = 0; i < bins; i++){
		unsigned short st[3];
		double r1 = i * step;
		for(j = 0; j <= bins; j++){
			double a = j / bins * M_PI;
			double s = sin(a/2.0), c = cos(a/2.0);

			for(k = 0; k < points_per_bin; k++){
				int bin, content_i, content_j, content_k;
				double x1[3], x2[3], x3[3], v1[3], v2[3];

				/* generate a point inside our volume */
				x1[0] = random_uniform_r(st, 0, mesh.size[0]);
				x1[1] = random_uniform_r(st, 0, mesh.size[1]);
				x1[2] = random_uniform_r(st, 0, mesh.size[2]);

				/* save contents of this location */
				content_i = (*mesh.voxel_content)((int)x1[0], (int)x1[1], (int)x1[2]) > 0;

				/* generate two random orthogonal vectors */
				random_unit_vector_r(st, v1);
				random_normal_r(st, v2, v1);

				x2[0] = x1[0] + r1 * (c * v1[0] + s * v2[0]);
				x2[1] = x1[1] + r1 * (c * v1[1] + s * v2[1]);
				x2[2] = x1[2] + r1 * (c * v1[2] + s * v2[2]);

				/* if x+v1 is outside our volume, discard it */
				if(x2[0] < 0 || x2[0] >= mesh.size[0] ||
				   x2[1] < 0 || x2[1] >= mesh.size[1] ||
				   x2[2] < 0 || x2[2] >= mesh.size[2] ){
					k = k-1;
					continue;
				}

				x3[0] = x1[0] + r1 * (c * v1[0] - s * v2[0]);
				x3[1] = x1[1] + r1 * (c * v1[1] - s * v2[1]);
				x3[2] = x1[2] + r1 * (c * v1[2] - s * v2[2]);

				/* if x+v2 is outside our volume, discard it */
				if(x3[0] < 0 || x3[0] >= mesh.size[0] ||
				   x3[1] < 0 || x3[1] >= mesh.size[1] ||
				   x3[2] < 0 || x3[2] >= mesh.size[2] ){
					k = k-1;
					continue;
				}

				/* save contents of other locations */
				content_j = (*mesh.voxel_content)((int)x2[0], (int)x2[1], (int)x2[2]) > 0;
				content_k = (*mesh.voxel_content)((int)x3[0], (int)x3[1], (int)x3[2]) > 0;

				bin = 2*content_i + 4*content_j + content_k;

#pragma omp atomic
				/* increase the right bin */
				Sijk[bin][i][j]++;
			}
		}
	}

	/* write result to output file */
	for(i = 0; i < bins; i++){
		for(j = 0; j <= bins; j++){
			fprintf(output, "%7.3f %7.3f", i*step*mesh.psize[0], (double)(180*j/bins));
			for(k = 0; k < 8; k++)
				fprintf(output," %f", (double) (Sijk[k][i][j])/points_per_bin);
			fprintf(output, "\n");
		}
		fprintf(output, "\n");
	}

	for(i = 0; i < 8; i++){
		for(j = 0; j < bins; j++)
			free(Sijk[i][j]);

		free(Sijk[i]);
	}
}
