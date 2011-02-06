#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mmf.h"

FILE *mmf_init(char *mmf_file,
               unsigned long *rows, unsigned long *cols, unsigned long *nz)
{
	char input[4096];
	int ret;
	FILE *f;

	if ( !(f = fopen(mmf_file, "r")) ){
		perror("mmf_init: fopen");
		exit(1);
	}

	while ( fgets(input, 4096, f)){
		if (input[0] == '#'){
			continue;
		}
		ret = sscanf(input, "%lu %lu %lu", rows, cols, nz);
		if ( ret != 3){
			fprintf(stderr, "mmf header error: sscanf");
			exit(1);
		}
		break;
	}

	return f;
}

int mmf_get_next(FILE *mmf, unsigned long *row, unsigned long *col, double *val)
{
	char input[4096];

	if ( !fgets(input, 4096, mmf) ){
		return 0;
	}
	if ( sscanf(input,  "%lu %lu %lf",  row, col, val) != 3 ){
		fprintf(stderr, "mmf file error: sscanf failed\n");
		fprintf(stderr, "**input:\n%s\n**\n", input);
		exit(1);
	}
	(*row)--; (*col)--;

	return 1;
}

int mmf_get_next_vstr(FILE *mmf, unsigned long *row, unsigned long *col, char **val)
{
	int ret;
	char input[4096], val_buf[4096];

	if ( !fgets(input, 4096, mmf) ){
		return 0;
	}
	ret = sscanf(input,  "%lu %lu %s",  row, col, val_buf);
	if (ret != 3){
		fprintf(stderr, "mmf file error: sscanf failed (%d) \n", ret);
		fprintf(stderr, "**input: %s**\n", input);
		exit(1);
	}
	(*row)--; (*col)--;

	int len = strlen(val_buf);
	val_buf[len] = '\n';
	val_buf[len+1] = '\0';
	*val = val_buf;

	return 1;
}
