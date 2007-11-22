/*
 * Program to query/enable/disable the hardware prefetcher
 * kkourt@cslab.ece.ntua.gr
 * References: 
 * IA32_MISC_ENABLE MSR in IA32 Vol. 3b
 */


#ifndef _LARGEFILE64_SOURCE
#define _LARGEFILE64_SOURCE
#endif

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <string.h>
#include <inttypes.h>

#define IA32_MISC_MSR 0x1a0

static inline uint64_t msr_set_bit(uint64_t val, int bit)
{
	return (val | (1<<bit));
}

static inline uint64_t msr_unset_bit(uint64_t val, int bit)
{
	uint64_t m = ~(1<<bit);
	return (val & m);
}

static inline unsigned msr_get_bit(uint64_t val, int bit)
{
	return ((val & (1<<bit)) != 0);
}

static struct {
	int       bit;
	char      *id;
	char      *desc;
} ia32_misc_bits[]  = {
	{ .bit = 7,  .id = "perfmon",    .desc = "Perfomance Monitoring Enabled"},
	{ .bit = 9,  .id = "pfetch-dlp", .desc = "Hardware Prefetch Disabled"},
	{ .bit = 19, .id = "pfetch-acl", .desc = "Adjacent Cache Line Prefetch Disabled"},
	{ 0 }
};

static void print_features()
{
	int i=0;
	
	printf("Features:\n");
	while ( ia32_misc_bits[i].id ){
		printf(" %-12s  %-40s\n", ia32_misc_bits[i].id,ia32_misc_bits[i].desc);
		i++;
	}

}

static int get_feature_bit(char *f)
{
	int i=0;
	char *id;

	while ( (id = ia32_misc_bits[i].id) ) {
		if ( strcmp(id, f) == 0){
			return (ia32_misc_bits[i].bit);
		}
		i++;
	}

	printf("error: cant find feature: %s\n", f);
	exit(1);
	return -1;
}

static void do_query(uint64_t val)
{
	char *id;
	char *dsc;
	int bit, i;

	for (i=0 ; (id =ia32_misc_bits[i].id)  ;i++){
		bit = ia32_misc_bits[i].bit;
		dsc = ia32_misc_bits[i].desc;

		printf("id:%-12s desc:%-40s  val: %d\n", id, dsc, msr_get_bit(val,bit));
	}
}


static void help(char *p)
{
	printf("Usage:\n");
	printf("%s query <cpu>\n", p);
	printf("%s (get|set|unset) <cpu> <feature>\n", p);
	print_features();
}

enum {
	OP_UNDF,
	OP_QUERY,
	OP_GET,
	OP_SET,
	OP_UNSET
};

static long parse_int(char *s)
{
	long ret;
	char *endptr;

	ret = strtol(s, &endptr, 10);
	if ( *endptr != '\0' ) {
		printf("parse error: '%s' is not a number\n", s);
		exit(1);
	}

	return ret;
}


int main(int argc, char **argv)
{
	int cpu=-1, msr_fd, ret, bit=-1;
	char *feat=NULL, msr_dev[64];
	unsigned op = OP_UNDF;
	uint64_t oldval, newval;

	if ( (argc == 3) && (strcmp(argv[1], "query") == 0) ){
		op  = OP_QUERY;
		cpu = parse_int(argv[2]);
	} else if ( (argc == 4) && (strcmp(argv[1], "get") == 0) ) {
		op   = OP_GET;
		cpu  = parse_int(argv[2]);
		feat = argv[3];
	} else if ( (argc == 4) && (strcmp(argv[1], "set") == 0) ) {
		op   = OP_SET;
		cpu  = parse_int(argv[2]);
		feat = argv[3];
	} else if ( (argc == 4) && (strcmp(argv[1], "unset") == 0) ) {
		op   = OP_UNSET;
		cpu  = parse_int(argv[2]);
		feat = argv[3];
	} 

	if ( op == OP_UNDF){
		help(argv[0]);
		exit(1);
	}
	
	snprintf(msr_dev,64,"/dev/cpu/%d/msr",cpu);

	msr_fd = open(msr_dev,O_RDWR);
	if ( msr_fd < 0){
		perror(msr_dev);
		exit(1);
	}

	ret = lseek(msr_fd, IA32_MISC_MSR, SEEK_SET);
	if ( ret == (off_t)-1){	
		perror("lseek");
		exit(1);
	}

	ret = read(msr_fd, &oldval, 8);
	if ( ret != 8 ) {
		perror("read");
		exit(1);
	}
		
	if ( feat )
		bit = get_feature_bit(feat);
	
	switch (op){
		case OP_QUERY:
		do_query(oldval);
		return 0;

		case OP_GET:
		printf("%d\n", msr_get_bit(oldval,bit));
		return 0;

		case OP_SET:
		newval = msr_set_bit(oldval,bit);
		break;

		case OP_UNSET:
		newval = msr_unset_bit(oldval,bit);
		break;
	}

	ret = lseek(msr_fd, IA32_MISC_MSR, SEEK_SET);
	if ( ret == (off_t)-1){	
		perror("lseek");
		exit(1);
	}

	ret = write(msr_fd, &newval, 8);
	if ( ret != 8){
		perror("write");
		exit(1);
	}

	return 0;
}
