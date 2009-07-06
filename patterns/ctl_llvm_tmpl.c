#include <stdio.h>
#include <inttypes.h>
#include <assert.h>

#include "ctl_ll.h"

typedef void llvm_jit_hook_t(void);

uint64_t u8_get(uint8_t **ctl)
{
	uint8_t ret = **ctl;
	(*ctl)++;

	return (uint64_t)ret;
}

uint64_t u16_get(uint8_t **ctl)
{
	uint16_t ret, **u16;

	u16 = (uint16_t **)ctl;
	ret = **u16;
	(*u16)++;

	return (uint64_t)ret;
}

uint64_t u32_get(uint8_t **ctl)
{
	uint32_t ret, **u32;

	u32 = (uint32_t **)ctl;
	ret = **u32;
	(*u32)++;

	return (uint64_t)ret;
}

uint64_t u64_get(uint8_t **ctl)
{
	uint64_t ret, **u64;

	u64 = (uint64_t **)ctl;
	ret = **u64;
	(*u64)++;

	return ret;
}

uint64_t ul_get(uint8_t **ctl)
{
	unsigned long ret;

	ret = u8_get(ctl);
	if (ret <= 127)
		goto end;

	unsigned shift = 7;
	unsigned long uc;
	ret -= 128;
	for (;;){
		uc = u8_get(ctl);
		if (uc <= 127){
			ret += (uc<<shift);
			break;
		}
		uc -= 128;
		ret += (uc<<shift);
		shift += 7;
	}
end:
	return ret;
}

#define ALIGN(buf,a) (void *) (((unsigned long) (buf) + (a-1)) & ~(a-1))
void align_ptr(uint8_t **ctl, int align)
{
	*ctl = ALIGN(*ctl, align);
}


void do_ctl_u8(uint8_t **ctl, uint64_t *x_indx, uint8_t size)
{
	int i;
	printf("x_indx:%lu\n", *x_indx);
	for (i=1; i<size; i++){
		*x_indx += u8_get(ctl);
		printf("x_indx:%lu\n", *x_indx);
	}
}

void do_ctl_u16(uint8_t **ctl, uint64_t *x_indx, uint8_t size)
{
	int i;
	for (i=1; i<size; i++){
		*x_indx += u16_get(ctl);
	}
}

llvm_jit_hook_t __new_row_hook __attribute__ ((annotate ("llvm_hook") ));
llvm_jit_hook_t __body_hook;
#if 0
void __attribute__ ((annotate ("llvm_hook") )) foo(void)
{
}
#endif

void
__llvm_annotate(void *val, char *descr, char *file, uint32_t line)
asm("llvm.var.annotation");

#define llvm_annotate(var, name) \
      __llvm_annotate(var, name, __FILE__, __LINE__)

void inline fail()
{
	assert(0);
}

void print_yx(uint64_t y, uint64_t x)
{
	printf("%lu %lu\n", y+1, x+1);
}

void ctl_decode_template(uint8_t *ctl, unsigned long ctl_size)
{
	uint8_t *ctl_end;
	uint64_t x_indx, y_indx;
	uint8_t size, flags;

	llvm_annotate(&x_indx, "vars::x_indx");
	llvm_annotate(&y_indx, "vars::y_indx");
	llvm_annotate(&flags, "vars::flags");
	llvm_annotate(&size, "vars::size");
	llvm_annotate(&ctl, "vars::ctl");

	ctl_end = ctl + ctl_size;
	x_indx = 0;
	y_indx = 0;
	do {
		flags = u8_get(&ctl);
		size = u8_get(&ctl);
		if (test_bit(&flags, CTL_NR_BIT)){
			__new_row_hook();
			//y_indx++;
			//x_indx = 0;
		}

		//printf("x_indx before jmp: %lu\n", x_indx);
		x_indx += ul_get(&ctl);
		//printf("x_indx after jmp: %lu\n", x_indx);
		__body_hook();

	} while (ctl < ctl_end);
}
