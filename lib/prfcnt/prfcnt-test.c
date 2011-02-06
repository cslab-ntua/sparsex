#include "prfcnt.h"

#include <stdlib.h>
#include <stdio.h>
#include <sched.h>
#include <linux/types.h> /* __u32, __u64 */
#include <sys/types.h>
#include <unistd.h>



#define ALIGN_BOUND 4095
#define ALIGN(buf)  (char *) (((unsigned long) (buf) + ALIGN_BOUND) & ~(ALIGN_BOUND))


int main()
{
	char *buffy1, *buffy1_orig, *buffy2, *buffy2_orig;
	int  i, err;
	int cpu=0;
	const unsigned int buffsize = (16*1024*1024);
	prfcnt_t prfcnt;
	
	cpu_set_t cpu_set;

	CPU_ZERO(&cpu_set);
	CPU_SET(cpu,&cpu_set);

	err = sched_setaffinity(getpid(), sizeof(cpu_set_t), &cpu_set);
	if (err){
		perror("sched_setaffinity");
		exit(1);
	}
	
	buffy1_orig = malloc(buffsize+ALIGN_BOUND);
	buffy2_orig = malloc(buffsize+ALIGN_BOUND);

	buffy1 = ALIGN(buffy1_orig);
	buffy2 = ALIGN(buffy2_orig);

	prfcnt_init(&prfcnt, cpu,PRFCNT_FL_T0|PRFCNT_FL_T1);

	prfcnt_start(&prfcnt);
	for (i=0; i < buffsize; i++){
		buffy1[i] = i;
	}
	for (i=0; i < buffsize; i++){
		buffy2[i] = buffy1[i] + 1;
	}
	prfcnt_pause(&prfcnt);
	
	prfcnt_report(&prfcnt);
	prfcnt_shut(&prfcnt);

	return 0;
}
