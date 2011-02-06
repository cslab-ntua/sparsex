#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#include "dynarray.h"
#include "vector.h"
#include "tsc.h"

#define BUFLEN 512

#define MAX(a,b) (a > b ? a : b)

int main(int argc, char **argv)
{
	char buf[BUFLEN], *s, *endptr;
	dynarray_t *idx_da;
	long vsize=0, size=0, i,j, loops;
	double v=0;
	tsc_t tsc;

	if (argc < 2){
		printf("Usage %s <loops>\n", argv[0]);
		exit(1);
	}
	loops = atol(argv[1]);
	if (!loops){
		printf("Usage %s <loops>\n", argv[0]);
		exit(1);
	}

	idx_da = dynarray_create(sizeof(uint32_t), 4096);

	for (;;){
		long idx;
		uint32_t *idx_ptr;
		s = fgets(buf,BUFLEN-1,stdin);
		if (s == NULL){
			break;
		}
		s[strlen(s) - 1] = '\0';

		idx = strtol(s, &endptr, 10);
		if (*endptr != '\0'){
			fprintf(stderr, "parsing error: %s not a number\n", s);
			exit(1);
		}

		idx_ptr = dynarray_alloc(idx_da);
		*idx_ptr = (uint32_t)idx;

		vsize = MAX(idx, vsize);
	}

	size = dynarray_size(idx_da);
	printf("using size:%ld vsize:%ld\n", size, vsize);

	vector_double_t *vd1 = vector_double_create(vsize);
	uint32_t *vis = dynarray_destroy(idx_da);
	vector_double_init(vd1, -0.666);

	vector_double_t *vd2 = vector_double_create(size);
	vector_double_init(vd2, -0.666);


	tsc_init(&tsc);
	tsc_start(&tsc);
	for (j=0;j<loops;j++){
		__asm__ __volatile__ (" # inner loop start\n");
		for (i=0;i<size;i++){
			v += vd1->elements[vis[i]] * vd2->elements[i];
		}
		__asm__ __volatile__ (" # inner loop end\n");
	}
	tsc_pause(&tsc);
	tsc_shut(&tsc);

	const double secs = tsc_getsecs(&tsc);
	const uint64_t ticks = tsc_getticks(&tsc);
	const uint64_t elems = 2*loops*size;
	printf("secs:%lf elems/sec:%lf (M) elems/ticks:%lf\n", 
	        secs, (((double)elems)/secs)/(1024.0*1024.0), ((double)elems/(double)ticks));

	printf("result: %lf\n", v);
	return 0;
}
