#include <sys/types.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/wait.h>

#include "ext_prog.h"


struct extp_s {
	int     p_r[2]; /* pipe for reading from stdout of external program */
	int     p_w[2]; /* pipe for writing to external program */
	int     pid;    /* pid of external program */
	char    *ext_prog; /* external program */
};

static void do_child(extp_t *extp)
{
	int err;

	close(extp->p_r[0]);
	close(extp->p_w[1]);

	err = dup2(extp->p_r[1], STDOUT_FILENO);
	if (err == -1){
		perror("dup2");
		exit(1);
	}
	close(extp->p_r[1]);

	err = dup2(extp->p_w[0], STDIN_FILENO);
	if (err == -1){
		perror("dup2");
		exit(1);
	}
	close(extp->p_w[0]);

	execl("/bin/sh", "sh", "-c", extp->ext_prog, NULL);

	perror("execl");
	exit(1);
}

extp_t *extp_open(char *ext_prog, FILE **in, FILE **out)
{
	extp_t *extp;
	int err;

	extp = malloc(sizeof(extp_t));
	if (!extp){
		perror("malloc");
		exit(1);
	}
	
	extp->ext_prog = ext_prog;

	err = pipe(extp->p_r);
	if (err){
		perror("pipe");
		exit(1);
	}

	err = pipe(extp->p_w);
	if (err){
		perror("pipe");
		exit(1);
	}

	extp->pid = fork();
	if (extp->pid == -1){
		perror("fork");
		exit(1);
	}

	if ( extp->pid == 0){
		do_child(extp);
	} 

	close(extp->p_r[1]);
	close(extp->p_w[0]);

	*in = fdopen(extp->p_r[0], "r");
	*out = fdopen(extp->p_w[1], "w");
	if (*in == NULL || *out == NULL){
		perror("fdopen");
		exit(1);
	}
	
	return extp;
}

int extp_close(extp_t *extp, FILE *in, FILE *out)
{
	int status, err;
	fclose(in);
	fclose(out);
	err = waitpid(extp->pid, &status, 0);
	if ( err == -1){
		perror("waitpid");
		exit(1);
	}
	
	return status;
}

#ifdef EXTP_MAIN_TEST
#define BUFLEN 1024
int main(int argc, char **argv)
{
	FILE *in, *out;
	extp_t *extp;
	char buff[BUFLEN], *s;
	int err;

	extp = extp_open("bc -l -q", &in, &out);

	for (;;){
		s = fgets(buff, BUFLEN-1, stdin);
		if (s == NULL){
			break;
		}

		err = fputs(s, out);
		if (err == EOF){
			fprintf(stderr, "cant write into pipe\n");
			exit(1);
		}
		//fflush(out);

		s = fgets(buff, BUFLEN-1, in);
		if (s == NULL){
			break;
		}

		printf("-->%s", s);
	}

	extp_close(extp, in, out);
	
	return 0;
}
#endif
