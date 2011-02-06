#include <stdlib.h>
#include <stdio.h>

#include "tsc.h"

int main(int argc, char **argv)
{
	if (argc < 2){
		fprintf(stderr, "Usage: %s <ticks>\n", argv[0]);
		exit(1);
	}

	tsc_spinticks(atoll(argv[1]));
	return 0;
}
