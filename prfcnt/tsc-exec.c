#include <unistd.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/wait.h>

#include "tsc.h"

/*
 * XXX: May not work as expected
 * see FIXME
 */
int main(int argc, char **argv)
{
	pid_t pid;
	int status;
	tsc_t timer;
	char **new_argv;

	if ( argc < 2){
		printf("Usage: %s <cmd> (args)\n", argv[0]);
		exit(1);
	}
	
	new_argv = &argv[1];

	if ( (pid = fork()) < 0){
		perror("fork");
		exit(1);
	}

	tsc_init(&timer);
	/*
	 * FIXME: Can we do this better ? 
	 */
	if (pid) {
		tsc_start(&timer);
		wait(&status);
		tsc_pause(&timer);
	} else {
		execv(argv[1],new_argv);
		perror("execv");
		exit(1);
	}

	tsc_report(&timer);
	tsc_shut(&timer);

	return 0;
}
