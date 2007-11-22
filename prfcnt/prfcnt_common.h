#ifndef __PRFCNT_COMMON_H__
#define __PRFCNT_COMMON_H__

/*
 * prfcnt_common.h
 *
 * Common functions (at least for x86 like insruction sets)
 * for accessing performance counters via :
 *  - /dev/msr linux device
 *  - rdpmc instruction 
 */

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#ifndef _LARGEFILE64_SOURCE
#define _LARGEFILE64_SOURCE
#endif

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <sys/types.h>
#include <unistd.h>
#include <sys/stat.h>
#include <fcntl.h>

#include <linux/types.h> /* __u32, __u64 */

static inline int prfcnt_wrmsr(int msr_fd, unsigned int addr, __u64 val)
{
	off_t ret;

	ret = lseek64(msr_fd, addr, SEEK_SET);
	if ( ret == (off_t)-1){
		perror("lseek");
		exit(1);
	}

	ret = write(msr_fd, &val, 8);
	if ( ret != 8 ){
		perror("write");
		exit(1);
	}

	return 0;
}

static inline int prfcnt_rdmsr(int msr_fd, unsigned int addr, __u64 *val)
{
	off_t ret;

	ret = lseek64(msr_fd, addr, SEEK_SET);
	if ( ret == (off_t)-1){
		perror("lseek");
		exit(1);
	}

	ret = read(msr_fd, val, 8);
	if ( ret != 8 ){
		perror("read");
		exit(1);
	}

	return 0;
}

static inline int __prfcnt_init(int cpu, int *msr_fd)
{
	int  fd;
	char msr_dev[32];
	snprintf(msr_dev, 32, "/dev/cpu/%d/msr", cpu);

	fd = open(msr_dev, O_RDWR);
	if (fd < 0){
		perror(msr_dev);
		exit(1);
	}

	*msr_fd = fd;

	return 0;
}

static inline void __prfcnt_shut(int msr_fd)
{
	close(msr_fd);
}

static inline __u64 prfcnt_rdpmc(__u32 cntr_nr)
{
	__u32 hi,low;
	__u64 ret;

	__asm__ __volatile__ ("rdpmc" : "=a"(low), "=d"(hi) : "c"(cntr_nr));

	ret = hi;
	ret <<= 32;
	ret |= low;

	return ret;
}

static inline __u32 prfcnt_rdpmc32(__u32 cntr_nr)
{
	__u32 ret;

	cntr_nr |= (1<<31);

	__asm__ __volatile__ ("rdpmc": "=a"(ret) : "c"(cntr_nr));

	return ret;
}

#endif
