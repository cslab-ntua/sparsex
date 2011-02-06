#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include "vector.h"
#include "tsc.h"

#define BUFLEN 512

int main(int argc, char **argv)
{
	long size, i,j, loops;
	double v=0;
	tsc_t tsc;

	if (argc < 3){
		printf("Usage %s <size> <loops>\n", argv[0]);
		exit(1);
	}
	size = atol(argv[1]);
	loops = atol(argv[2]);
	if (!size || !loops){
		printf("Usage %s <size> <loops>", argv[0]);
		exit(1);
	}

	vector_double_t *vd1 = vector_double_create(size);
	vector_double_t *vd2 = vector_double_create(size);
	vector_double_init(vd1, -0.666);
	vector_double_init(vd2, -0.666);

	tsc_init(&tsc);
	tsc_start(&tsc);
	for (j=0;j<loops;j++){
		__asm__ __volatile__ (" # inner loop start\n");
		for (i=0;i<size;i++){
			v += vd1->elements[i] * vd2->elements[i];;
		}
		__asm__ __volatile__ (" # inner loop end\n");
	}
	tsc_pause(&tsc);
	tsc_shut(&tsc);

	const double secs = tsc_getsecs(&tsc);
	const uint64_t ticks = tsc_getticks(&tsc);
	const uint64_t elems = loops*size*2;
	printf("secs:%lf elems/sec:%lf (M) elems/ticks:%lf\n", 
	        secs, (((double)elems)/secs)/(1024.0*1024.0), ((double)elems/(double)ticks));

	printf("result: %lf\n", v);
	return 0;
}
