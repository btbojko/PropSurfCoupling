/* ROCSTAT - slice.c */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "amira.h"

extern FILE* output;

void make_slice(char* slice){
	register int i, j, k;
	char *stype, *snum;
	int   slice_number;

	stype = strtok(slice,":\n");
	snum  = strtok((char *)NULL,":\n");

	/* some error checking */
	if(stype == NULL || snum == NULL){
		fprintf(stderr, "Error: %s:%s is not a valid argument to option -s.\n", stype, snum);
		exit(1);
	} else {
		slice_number = strtol(snum, (char **)NULL, 10);
	}

	if(strcmp(stype, "xz") == 0){
		j = slice_number - 1;

		if(j < 0 || j > mesh.size[1]){
			fprintf(stderr,"Slice number out of range: %d (should be 1 to %lu)\n", slice_number, mesh.size[1]);
			exit(1);
		}

		for(k = 0; k < mesh.size[2]; k++){
			for(i = 0; i < mesh.size[0]; i++)
				fprintf(output,"%hu ", (*mesh.voxel_content)(i,j,k));
			fprintf(output,"\n");
		}

	} else if(strcmp(stype,"xy") == 0){
		k = slice_number - 1;

		if(k < 0 || k > mesh.size[2]){
			fprintf(stderr,"Slice number out of range: %d (should be 1 to %lu)\n",slice_number, mesh.size[2]);
			exit(1);
		}

		for(j = 0; j < mesh.size[1]; j++){
			for(i = 0; i < mesh.size[0]; i++)
				fprintf(output,"%hu ", (*mesh.voxel_content)(i,j,k));
			fprintf(output,"\n");
		}

	} else if(strcmp(stype,"yz") == 0){
		i = slice_number - 1;

		if(i < 0 || i > mesh.size[0]){
			fprintf(stderr,"Slice number out of range: %d (should be 1 to %lu)\n", slice_number, mesh.size[0]);
			exit(1);
		}

		for(k = 0; k < mesh.size[2]; k++){
			for(j = 0; j < mesh.size[1]; j++)
				fprintf(output,"%hu ", (*mesh.voxel_content)(i,j,k));
			fprintf(output,"\n");
		}

	} else {
		fprintf(stderr,"Slice type \"%s\" not recognized.\n", stype);
		exit(1);
	}
}
