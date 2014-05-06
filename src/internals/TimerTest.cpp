#include "timer.hpp"
#include <stdio.h>

int main()
{
	timing::Timer timer;

	printf("sleeping for 1 second and reporting:\n");
	timer.Start();
	sleep(1);
	timer.Pause();
	printf("timed %lf secs\n", timer.ElapsedTime());

    timer.Clear();
	printf("sleeping 2 times for 1 second and reporting:\n");
	timer.Start();
	sleep(1);
	timer.Pause();

	sleep(1);

	timer.Start();
	sleep(1);
	timer.Pause();
	printf("timed %lf secs\n", timer.ElapsedTime());
	
	return 0;
}
