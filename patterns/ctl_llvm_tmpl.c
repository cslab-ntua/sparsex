#include <stdio.h>
#include <inttypes.h>
#include <assert.h>

#include "ctl_ll.h"

typedef void llvm_jit_hook_t(void);

uint8_t u8_get(uint8_t **ctl)
{
	uint8_t ret = **ctl;
	(*ctl)++;

	return ret;
}

uint16_t u16_get(uint8_t **ctl)
{
	uint16_t ret, **u16;

	u16 = (uint16_t **)ctl;
	ret = **u16;
	(*u16)++;

	return ret;
}

uint32_t u32_get(uint8_t **ctl)
{
	uint32_t ret, **u32;

	u32 = (uint32_t **)ctl;
	ret = **u32;
	(*u32)++;

	return ret;
}

uint64_t u64_get(uint8_t **ctl)
{
	uint64_t ret, **u64;

	u64 = (uint64_t **)ctl;
	ret = **u64;
	(*u64)++;

	return ret;
}

unsigned long ul_get(uint8_t **ctl)
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

llvm_jit_hook_t __new_row_hook __attribute__ ((annotate ("llvm_hook") ));
llvm_jit_hook_t __body_hook;

void __attribute__ ((annotate ("llvm_hook") )) foo(void)
{
}


void
__llvm_annotate(void *val, char *descr, char *file, uint32_t line)
asm("llvm.var.annotation");

void ctl_decode_template(uint8_t *ctl, unsigned long ctl_size)
{
	uint8_t *ctl_end;
	uint64_t x_indx, y_indx;

	__llvm_annotate(&x_indx, "vars::x_indx", __FILE__, __LINE__);
	__llvm_annotate(&y_indx, "vars::y_indx", __FILE__, __LINE__);

	ctl_end = ctl + ctl_size;
	x_indx = 0;
	y_indx = 0;
	do {
		if (test_bit(ctl, CTL_NR_BIT)){
			__new_row_hook();
		}

		x_indx += ul_get(&ctl);

		__body_hook();

	} while (ctl < ctl_end);
}
