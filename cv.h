#ifndef __VC_H__
#define __VC_H__

#include <stdlib.h>
#include <stdio.h>

#define DECLARE_UC_GET_INT(uc_pointer) \
	static  int uc_get_int(void)\
	{\
		static unsigned long uc_idx=0;\
		int val;\
		unsigned int tmp, neg;\
	\
		__asm__ __volatile__ (\
			"xorl %[val], %[val] \n\t"\
			"movzbl (%[uc_ptr], %[uc_idx]),%[tmp]\n\t"\
			"incq %[uc_idx]\n\t"\
			"btrl $7, %[tmp] \n\t"\
			"sbbl %[neg], %[neg] \n\t"\
			"btrl $6, %[tmp] \n\t"\
			"jnc 1f\n\t" \
	\
			"orl %[tmp], %[val]\n\t"\
			"movzbl (%[uc_ptr], %[uc_idx]),%[tmp]\n\t"\
			"incq %[uc_idx]\n\t"\
			"shl $6, %[tmp] \n\t"\
			"btrl $13, %[tmp] \n\t"\
			"jnc 1f\n\t" \
	\
			"orl %[tmp], %[val]\n\t"\
			"movzbl (%[uc_ptr], %[uc_idx]),%[tmp]\n\t"\
			"incq %[uc_idx]\n\t"\
			"shl $13, %[tmp] \n\t"\
			"btrl $20, %[tmp] \n\t"\
			"jnc 1f\n\t" \
	\
			"orl %[tmp], %[val]\n\t"\
			"movzbl (%[uc_ptr], %[uc_idx]),%[tmp]\n\t"\
			"incq %[uc_idx]\n\t"\
			"shl $20, %[tmp] \n\t"\
			"btrl $27, %[tmp] \n\t"\
	\
			"1: \n\t"\
			"orl %[tmp], %[val]\n\t"\
			"movl %[val], %[tmp]\n\t"\
			"negl %[tmp] \n\t"\
			"testl %[neg], %[neg] \n\t"\
			"cmovne %[tmp], %[val] \n\t"\
			: [uc_idx] "+r" (uc_idx),\
			  [val] "=r" (val),\
			  [tmp] "=&r" (tmp),\
			  [neg] "=&r" (neg)\
			: [uc_ptr] "r" (uc_pointer)\
		);\
	\
		return val;\
	}

#define DECLARE_UC_PUT_INT \
	static void uc_put_int(long i, unsigned char **uc_pptr)\
	{\
		unsigned char *uc_ptr = *uc_pptr;\
		long _i = i;\
		if (i < 0){\
			*uc_ptr = (1<<7);\
			i = -i;\
		}\
		*uc_ptr++ |= ( i & ((1<<6) - 1) );\
		i >>= 6;\
		if ( !i ) return;\
		*(uc_ptr-1) |= (1<<6);\
	\
		*uc_ptr++ = ( i & ((1<<7) - 1) );\
		i >>= 7;\
		if ( !i ) return;\
		*(uc_ptr-1) |= (1<<7);\
	\
		*uc_ptr++ = ( i & ((1<<7) - 1) );\
		i >>= 7;\
		if ( !i ) return;\
		*(uc_ptr-1) |= (1<<7);\
	\
		*uc_ptr++ = ( i & ((1<<7) - 1) );\
		i >>= 7;\
		if ( !i ) return;\
		*(uc_ptr-1) |= (1<<7);\
		\
		printf("%ld does not fit\n", _i);\
		exit(1);\
	}

#define MAX_INT_SUPPORTED 134217727

#ifdef VC_TEST
#include <inttypes.h>

#include "tsc.h"

int main(int argc, char **argv)
{
	#if 0
	int val, new_val;
	uc_ptr = foo;
	
	unsigned char foo[4] = {0, 0, 0, 0};
	if (argc < 2){
		printf("usage: %s <int>\n", argv[0]);
		exit(1);
	}

	val = strtoimax(argv[1], &error, 10);
	if ( *error != '\0' ){
		printf("error: %s is not an int\n", argv[1]);
		exit(1);
	}
	uc_put_int(val);

	int i;
	for (i=0; i<4; i++){
		printf("cu-- %d, %lu\n", i, (unsigned long)foo[i]);
	}

	uc_ptr = foo;
	new_val = uc_get_int();
	printf("%d %d\n", val, new_val);
	#endif
	unsigned char *foo = malloc(MAX_INT_SUPPORTED*4);
	long i = 0;
	if ( !foo ){
		perror("malloc failed");
		exit(1);
	}
	uc_ptr = foo;

	for (i=0; i<MAX_INT_SUPPORTED; i++){
		uc_put_int((int)(i%64));
	}
	
	uc_ptr = foo;
	tsc_t tsc;
	unsigned long ret=0;
	tsc_init(&tsc);
	tsc_start(&tsc);
	for (i=0; i<MAX_INT_SUPPORTED; i++){
		ret += uc_get_int();
	}
	tsc_pause(&tsc);
	tsc_report(&tsc);

	printf("ret=%lu\n", ret);


	return 0;
}

#endif /* VC_TEST */

#endif /* __VC_H__ */
