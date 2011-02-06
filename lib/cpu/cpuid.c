/*
 * find out cpuid info
 * kkourt@cslab.ece.ntua.gr
 */

#define _LARGEFILE64_SOURCE

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <string.h>

#include <linux/types.h>

typedef void cpuid_fn_t(unsigned int index, unsigned int *results);

static int fd;
cpuid_fn_t *cpuid_fn;

void usage()
{
	printf("Usage: cpuid <cpu (def=0)> <cmd (def=print)>\n");
	printf("Commands:\n");
	printf("\tprint         : read index from stdin and print results\n");
	printf("\tbrand_string  : Print processor brand string\n");
	printf("\tarch_perf     : Architectural Performance Model\n");
}

#define bitfield_read(data, start, end) ((data>>start) & (~(~0U<<(end-start+1))))

void cpuid_device(unsigned int index, unsigned int *results)
{
	off64_t       ret, seek;

	seek = (off64_t)index;
	ret = lseek(fd, seek, SEEK_SET);
	if (ret == (off64_t)-1){
		perror("lseek");
		return;
	}

	ret = read(fd, results, 16);
	if (ret == -1){
		perror("read");
		return;
	}
}

void cpuid_asm(unsigned int index, unsigned int *results)
{
	unsigned int eax, ebx, ecx, edx;

	#if defined(__i386__) || defined(__x86_64__)
	__asm__ __volatile("cpuid"
	                    : "=a"(eax), "=b"(ebx), "=c"(ecx), "=d"(edx)
			    : "a" (index));
	#else
	eax = ebx = ecx = edx = 0;
	#endif

	results[0] = eax;
	results[1] = ebx;
	results[2] = ecx;
	results[3] = edx;
}

void cpuid_print()
{
	char buff[64];
	unsigned int results[4];
	int index;
	
	//printf("please type input: ");
	while ( fgets(buff,63,stdin) != NULL) {
		
		index = strtoul(buff,NULL,0);
		cpuid_fn(index, results);
		printf("cpuid\t (0x%x) eax:0x%x ebx:0x%x ecx:0x%x edx:0x%x\n",
		       index, results[0], results[1], results[2], results[3]);
	}
}

/*
 * IA32-2a CPUID instruction 
 */
void cpuid_brand_string()
{
	unsigned int results[4];
	char str[48];
	unsigned int max;
	
	cpuid_fn(0x80000000, results);

	max = results[0];
	if ( !( (max & 0x80000000) && (max > 0x80000004) ) ){
		printf("Processor Brand String Not Supported\n");
		return;
	}
	
	cpuid_fn(0x80000002, (void *)str);
	cpuid_fn(0x80000003, (void *)(str+16));
	cpuid_fn(0x80000004, (void *)(str+32));

	printf("%s\n", str);
}

void cpuid_arch_perf()
{
	unsigned int results[4];
	unsigned int eax, ebx, v;

	cpuid_fn(0x0a, results);

	eax = results[0];
	printf("eax=0x%x\n", eax);
	printf("Version ID                      %3u\n", v = bitfield_read(eax,0 ,7 ));
	if (!v) {
		printf("This processor doesn't support architectural perf. cntrs\n");
		return;
	}
	printf("Perf. mon. cntr / logical cpu   %3u\n", bitfield_read(eax,8 ,15)); 
	printf("Perf. mon. cntr Bit width       %3u\n", bitfield_read(eax,16,23));
	printf("Length of EBX bit vector        %3u\n", bitfield_read(eax,24,31));

	printf("\n");
	ebx = results[1];
	printf("Core cycle             %u\n", !bitfield_read(ebx,0,0));
	printf("Instr retired          %u\n", !bitfield_read(ebx,1,1));
	printf("Ref. cycles            %u\n", !bitfield_read(ebx,2,2));
	printf("Last Level Cache Ref.  %u\n", !bitfield_read(ebx,3,3));
	printf("Last Level Cache Miss  %u\n", !bitfield_read(ebx,4,4));
	printf("Branch Instr. retired  %u\n", !bitfield_read(ebx,5,5));
	printf("Branch Miss retired    %u\n", !bitfield_read(ebx,6,6));
}

int main(int argc, char **argv)
{
	char cpuname[32];
	char *cmd = NULL;
	int cpu=0;

	while ( argc-- > 1 ){
		if ( argv[1][0] >= '0' && argv[1][0] <= '9' ){
			cpu = atoi(argv[1]);
		} else {
			cmd = argv[1];
		}

		argv++;
	}

	snprintf(cpuname,20,"/dev/cpu/%d/cpuid", cpu);
	
	fd = open(cpuname, O_RDONLY);
	if (fd == -1){
		fprintf(stderr, "could not open %s \n", cpuname);
		cpuid_fn = cpuid_asm;
	} else {
		cpuid_fn = cpuid_device;
	}
	
	if ( !cmd || !strcmp(cmd, "print") ) {
		cpuid_print();
	} else if ( !strcmp(cmd, "brand_string") ) {
		cpuid_brand_string();
	} else if ( !strcmp(cmd, "arch_perf") ) {
		cpuid_arch_perf();
	} else {
		usage();
	}

	if ( fd != -1) 
		close(fd);

	return 0;
}
