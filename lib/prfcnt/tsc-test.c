#include <stdio.h>
#include <unistd.h>

#include "tsc.h"

int main()
{
	tsc_t timer;

	tsc_init(&timer);

	printf("sleeping for 1 second and reporting:\n");
	tsc_start(&timer);
	sleep(1);
	tsc_pause(&timer);
	tsc_report(&timer);

	printf("sleeping 2 times for 1 second and reporting:\n");
	tsc_init(&timer);

	tsc_start(&timer);
	sleep(1);
	tsc_pause(&timer);

	sleep(1);

	tsc_start(&timer);
	sleep(1);
	tsc_pause(&timer);
	
	tsc_report(&timer);
	tsc_shut(&timer);
	
	return 0;
}
