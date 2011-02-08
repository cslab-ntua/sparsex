#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <inttypes.h> /* PRIu64 */

#include "mmf.h"

FILE *mmf_init(char *mmf_file,
               uint64_t *rows, uint64_t *cols, uint64_t *nz)
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
		ret = sscanf(input, "%" PRIu64 " %" PRIu64 " %" PRIu64, rows, cols, nz);
		if ( ret != 3){
			fprintf(stderr, "mmf header error: sscanf");
			exit(1);
		}
		break;
	}

	return f;
}

int mmf_get_next(FILE *mmf, uint64_t *row, uint64_t *col, double *val)
{
	char input[4096];

	if ( !fgets(input, 4096, mmf) ){
		return 0;
	}
	if ( sscanf(input,  "%"PRIu64" %"PRIu64" %lf",  row, col, val) != 3 ){
		fprintf(stderr, "mmf file error: sscanf failed\n");
		fprintf(stderr, "**input:\n%s\n**\n", input);
		exit(1);
	}
	(*row)--; (*col)--;

	return 1;
}

int mmf_get_next_vstr(FILE *mmf, uint64_t *row, uint64_t *col, char **val)
{
	int ret;
	char input[4096], val_buf[4096];

	if ( !fgets(input, 4096, mmf) ){
		return 0;
	}
	ret = sscanf(input,  "%"PRIu64" %"PRIu64" %s",  row, col, val_buf);
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
