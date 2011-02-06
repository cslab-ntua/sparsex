#include <sched.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libgen.h>
#include <time.h>

#include "tsc.h"
#include "prfcnt.h"
#include "matrix.h"

#define LOOPS 128
#define MAT_SIZE 1024


int pollute_cache()
{
	char *b;
	int  i, sum;
	const int size = (CACHE_BYTES << 1);

	b = malloc(size);
	if (!b){
		perror("malloc");
		exit(1);
	}

	for (i=0 ; i < size; i++)
		b[i] = i;

	for (sum=0, i=0 ; i < size; i++)
		sum += b[i];

	free(b);

	return sum;
}

static void do_setaffinity(cpu_set_t *mask)
{
	int err;

	err = sched_setaffinity(0, sizeof(cpu_set_t), mask);
	if (err){
		perror("sched_setaffinity");
		exit(1);
	}
}

static void vxv_mul(vector_t *a, elem_t *b, elem_t *out)
{
	unsigned long i;

	elem_t *v1 = a->elements;
	elem_t v2 = *b;
	unsigned long size = a->size;

	for (i=0; i<size; i++){
		*out += (v1[i]*v2);
	}

}

int main(int argc, char **argv)
{
	int fd, i, ret;
	int cpu=0;
	vector_t *v;
	elem_t y, x;
	tsc_t tsc;
	prfcnt_t prfcnt;
	cpu_set_t cpu_mask;
	unsigned long size;

	size = (argc<2) ? MAT_SIZE: atoi(argv[1]);

	CPU_ZERO(&cpu_mask);
	CPU_SET(cpu, &cpu_mask);
	do_setaffinity(&cpu_mask);

	v = vector_create(size);
	y = x = 0;

	srand(time(NULL));
	tsc_init(&tsc);
	prfcnt_init(&prfcnt,cpu,PRFCNT_FL_T0|PRFCNT_FL_T1);
	for (i=0; i<LOOPS ;  i++){
		vector_init_rand_range(v, (elem_t)-1000, (elem_t)1000);
		pollute_cache();
		prfcnt_start(&prfcnt);
		tsc_start(&tsc);
		vxv_mul(v, &x, &y);
		tsc_pause(&tsc);
		prfcnt_pause(&prfcnt);
	}
	//tsc_report(&tsc);
	tsc_shut(&tsc);
	float mflops = (((float)(size*2*LOOPS))/tsc_getsecs(&tsc))/(1000*1000);
	printf("%lu %f\n", size, mflops);
	prfcnt_report(&prfcnt);

	return 0;
}


