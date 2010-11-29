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

void print_yxv(uint64_t y, uint64_t x, double v)
{
	printf("%lu %lu %le\n", y+1, x+1, v);
}

#if 0
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
#endif

#define ELEM_TYPE double
#include "vector.h"
#include "csx.h"

void csx_spmv_template(void *spm, vector_double_t *in, vector_double_t *out)
{
    csx_double_t *csx = (csx_double_t *)spm;
    double *x;
    double *y;
    double *v = csx->values;
    double *myx;
    double yr = 0;
    uint8_t *ctl = csx->ctl;
    uint8_t *ctl_end = ctl + csx->ctl_size;
    uint64_t y_indx = csx->row_start;
    uint8_t size, flags;

    //printf("csx->ctl: %p\n", csx->ctl);

    x = in  ? in->elements : NULL;
    y = out ? out->elements: NULL;
    myx = x;

    llvm_annotate(&yr, "spmv::yr");
    llvm_annotate(&myx, "spmv::myx");
    llvm_annotate(&x, "spmv::x");
    llvm_annotate(&y, "spmv::y");
    llvm_annotate(&y_indx, "spmv::y_indx");
    llvm_annotate(&v, "spmv::v");
    llvm_annotate(&ctl, "spmv::ctl");
    llvm_annotate(&size, "spmv::size");
    llvm_annotate(&flags, "spmv::flags");
    do {
        //printf("ctl:%p\n", ctl);
        flags = *ctl++;
        size = *ctl++;
        //printf("size=%d\n", size);
        if (test_bit(&flags, CTL_NR_BIT)){
            __new_row_hook();
            myx = x;
            yr = 0;
            //y[y_indx] = yr;
        }
        //printf("x_indx before jmp: %lu\n", myx - x);
        myx += ul_get(&ctl);
        //printf("x_indx after jmp: %lu\n", myx - x);
        __body_hook();
        //printf("x_indx at end: %lu\n", myx - x);
    } while (ctl < ctl_end);
    y[y_indx] += yr;
}
