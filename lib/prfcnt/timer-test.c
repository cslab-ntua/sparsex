#include <stdio.h>
#include <unistd.h>

#include "timer.h"

int main()
{
	xtimer_t timer;

	timer_init(&timer);

	printf("sleeping for 1 second and reporting:\n");
	timer_start(&timer);
	sleep(1);
	timer_pause(&timer);
	printf("timed %lf secs\n", timer_secs(&timer));

	printf("sleeping 2 times for 1 second and reporting:\n");
	timer_init(&timer);

	timer_start(&timer);
	sleep(1);
	timer_pause(&timer);

	sleep(1);

	timer_start(&timer);
	sleep(1);
	timer_pause(&timer);
	printf("timed %lf secs\n", timer_secs(&timer));
	
	return 0;
}
