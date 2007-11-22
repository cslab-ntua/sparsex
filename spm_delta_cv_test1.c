#include <stdio.h>
#include <stdlib.h>

#include "spm_delta_cv.h"

int main(int argc, char **argv)
{
	char buf[16];
	long i1=0, i2=0;
	long c1, c2;
	int size;

	if ( argc > 1){
		i1 = atol(argv[1]);
		i2 = -i1;
	}

	for (;;){

		printf("checking %ld\n", i1);
		size = spm_delta_visize(i1);
		si_set(buf, i1, size);
		c1 = si_get(buf, size);
		if (i1 != c1){
			printf("set:%ld got:%ld (size:%d) \n", i1, c1, size);
			exit(1);
		}

		printf("checking %ld\n", i2);
		size = spm_delta_visize(i2);
		si_set(buf, i2, size);
		c2 = si_get(buf, size);
		if (i2 != c2){
			printf("set:%ld got:%ld (size:%d)\n", i2, c2, size);
			exit(1);
		}

		i1++; i2--;
	}
}
