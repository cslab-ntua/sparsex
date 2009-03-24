#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>

int main(int argc, char **argv)
{
	int vi;

	if (argc < 2){
		fprintf(stderr, "Usage: %s <range>\n", argv[0]);
		exit(1);
	}

	unsigned long distinct_vals = atol(argv[1]);
	if (!distinct_vals){
		fprintf(stderr, "invalid range: %s\n", argv[1]);
		exit(1);
	}
	if ( distinct_vals < 1UL<<(sizeof(uint8_t)*8) ){
		vi = 8;
	} else if ( distinct_vals < 1UL<<(sizeof(uint16_t)*8) ){
		vi = 16;
	} else if ( distinct_vals < 1UL<<(sizeof(uint32_t)*8) ){
		vi = 32;
	} else {
		fprintf(stderr, "range: %lu too big\n", distinct_vals);
		exit(1);
	}

	printf("%d\n", vi);
	return 0;
}
