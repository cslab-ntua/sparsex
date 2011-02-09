#include "prfcnt.h"
#include "tsc.h"

#include <unistd.h>
#include <stdlib.h>
#include <sched.h>
#include <sys/types.h>
#include <sys/wait.h>


/*
 * XXX: May not work as expected
 * see FIXME
 */

int main(int argc, char **argv)
{
	pid_t pid;
	int status;
	char **new_argv;
	prfcnt_t prfcnt;
	tsc_t timer;
	cpu_set_t cpu_set;
	int err;

	if ( argc < 2){
		printf("Usage: %s <cmd> (args)\n", argv[0]);
		exit(1);
	}

	new_argv = &argv[1];

	/*
	 * CPU affinity is inherited across a fork()
	 */
	CPU_ZERO(&cpu_set);
	CPU_SET(0,&cpu_set);

	err = sched_setaffinity(getpid(), sizeof(cpu_set_t), &cpu_set);
	if (err){
		perror("sched_setaffinity");
		exit(1);
	}

	if ( (pid = fork()) < 0){
		perror("fork");
		exit(1);
	}

	tsc_init(&timer);
	prfcnt_init(&prfcnt,0,PRFCNT_FL_T0|PRFCNT_FL_T1);

	/*
	 * FIXME: Is this efficient enough ? Could it be done better ?
	 */
	if (pid) {
		prfcnt_start(&prfcnt);
		tsc_start(&timer);
		wait(&status);
		tsc_pause(&timer);
		prfcnt_pause(&prfcnt);

	} else {

		execv(argv[1],new_argv);
		perror("execv");
		exit(1);
	}

	tsc_report(&timer);
	prfcnt_report(&prfcnt);

	return 0;
}
