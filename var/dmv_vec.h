#ifndef __DMV_VEC_H__
#define __DMV_VEC_H__

typedef double dv_t __attribute__ ((vector_size (16)));
typedef float  fv_t __attribute__ ((vector_size (16)));

typedef union dvec {
	dv_t dv;
	struct { double _d0, _d1; } d;
} dvec_t;

typedef union fvec {
	fv_t fv;
	struct { float _f0, _f1, _f2, _f3; } f;
} fvec_t;

#include "matrix.h"

#if ELEM_SIZE == 8
#define DMV_VEC_SIZE 2
static inline double dmv_vec(double *v, double *x, long nr)
{
	#if 0
	if ( ((unsigned long)v & 15) ){
		printf("v unaligned: %p\n", v);
		exit(1);
	}
	if ( ((unsigned long)x & 15) ){
		printf("x unaligned: %p\n", x);
		exit(1);
	}
	#endif

	unsigned long _cnt;
	register dvec_t vr, xr, yr;
	dv_t *v_vec = (dv_t *)v;
	dv_t *x_vec = (dv_t *)x;
	yr.d._d0 = yr.d._d1 = (double)0;
	__asm__ __volatile__(" # dmv_vec start");
	for ( _cnt=0; _cnt<nr; _cnt+=2){
		vr.dv = *v_vec++;
		xr.dv = *x_vec++;
		//printf("x0:   %-10.10lf x1:   %-10.10lf\n", xr.d._d0, xr.d._d1);
		//printf("v0:   %-10.10lf v1:   %-10.10lf\n", vr.d._d0, vr.d._d1);
		xr.dv *= vr.dv;
		//printf("mul0: %-10.10lf mul1: %-10.10lf\n", xr.d._d0, xr.d._d1);
		yr.dv += xr.dv;
		//printf("add0: %-10.10lf add1: %-10.10lf\n\n", yr.d._d0, yr.d._d1);
	}
	__asm__ __volatile__(" # dmv_vec end");
	return (yr.d._d0 + yr.d._d1);
} 

static inline double dmv_vec_u(double *v, double *x, long nr)
{
	register unsigned long cnt, not_aligned;

	register double yr_na = (double) 0;
	not_aligned = nr & (DMV_VEC_SIZE-1);
	for (cnt=0 ; cnt < not_aligned; cnt++){
		yr_na += v[cnt]*x[cnt];
	}
	//printf("na:%lu yr_na:%lf\n", not_aligned, yr_na);

	register dvec_t xreg, yreg, vreg;
	yreg.d._d0 = yreg.d._d1 = (double)0;
	__asm__ __volatile__ (" # start dmv_vec2");
	for ( ; cnt < nr ; cnt += DMV_VEC_SIZE){
		xreg.dv = __builtin_ia32_loadupd(x+cnt);
		//printf("x0:   %-10.10lf x1:   %-10.10lf\n", xreg.d._d0, xreg.d._d1);
		vreg.dv = __builtin_ia32_loadupd(v+cnt);
		//printf("v0:   %-10.10lf v1:   %-10.10lf\n", vreg.d._d0, vreg.d._d1);
		xreg.dv = __builtin_ia32_mulpd(vreg.dv, xreg.dv);
		//printf("mul0: %-10.10lf mul1: %-10.10lf\n", xreg.d._d0, xreg.d._d1);
		yreg.dv = __builtin_ia32_addpd(xreg.dv, yreg.dv);
		//printf("add0: %-10.10lf add1: %-10.10lf\n\n", yreg.d._d0, yreg.d._d1);
	}

	return (yreg.d._d0 + yreg.d._d1 + yr_na);
}

static inline double dmv_vec_va(double *v, double *x, long nr)
{
	register unsigned long cnt, not_aligned;
	register double yr_na = (double) 0;
	not_aligned = nr & (DMV_VEC_SIZE-1);
	for (cnt=0 ; cnt < not_aligned; cnt++){
		yr_na += v[cnt]*x[cnt];
	}

	register dvec_t xreg, yreg, vreg;
	dv_t *v_vec = (dv_t *)v;
	__asm__ __volatile__ (" # start dmv_vec2");
	for ( ; cnt < nr ; cnt += DMV_VEC_SIZE){
		xreg.dv = __builtin_ia32_loadupd(x+cnt);
		xreg.dv = __builtin_ia32_mulpd(*v_vec++, xreg.dv);
		yreg.dv = __builtin_ia32_addpd(xreg.dv, yreg.dv);
	}

	return (yreg.d._d0 + yreg.d._d1 + yr_na);
}

/* XXX: let the compiler do it */
#if 0
static inline double dmv_vec_ur(double *v, double *x, long nr)
{
	unsigned long _cnt;
	register dvec_t xr, yr;
	register dvec_t vr0, vr1, vr2, vr3, vr4, vr5, vr6, vr7;
	dv_t *v_vec = (dv_t *)v;
	dv_t *x_vec = (dv_t *)x;
	yr.d._d0 = yr.d._d1 = (double)0;
	__asm__ __volatile__(" # dmv_vec_ur start");
	for ( _cnt=0; _cnt<nr; _cnt+=2*8){

		vr0.dv = *(v_vec + 0);
		vr1.dv = *(v_vec + 1);
		vr2.dv = *(v_vec + 2);
		vr3.dv = *(v_vec + 3);
		vr4.dv = *(v_vec + 4);
		vr5.dv = *(v_vec + 5);
		vr6.dv = *(v_vec + 6);
		vr7.dv = *(v_vec + 7);
		v_vec += 8;

		vr0.dv *= *x_vec++;
		vr1.dv *= *x_vec++;
		vr2.dv *= *x_vec++;
		vr3.dv *= *x_vec++;
		vr4.dv *= *x_vec++;
		vr5.dv *= *x_vec++;
		vr6.dv *= *x_vec++;
		vr7.dv *= *x_vec++;

		vr0.dv += vr1.dv;
		vr2.dv += vr3.dv;
		vr4.dv += vr5.dv;
		vr6.dv += vr7.dv;

		vr0.dv += vr2.dv;
		vr4.dv += vr6.dv;

		yr.dv = vr0.dv + vr4.dv;
		
	}
	__asm__ __volatile__(" # dmv_vec_ur end");
	return (yr.d._d0 + yr.d._d1);
} 
#endif

#elif ELEM_SIZE == 4
#define DMV_VEC_SIZE 4
static inline float dmv_vec(float *v, float *x, long nr)
{
	unsigned long _cnt;
	register fvec_t vr, xr, yr;
	fv_t *v_vec = (fv_t *)v;
	fv_t *x_vec = (fv_t *)x;
	yr.f._f0 = yr.f._f1 = yr.f._f2 = yr.f._f3 = (float)0;
	__asm__ __volatile__(" # dmv_vec start");
	for ( _cnt=0; _cnt<nr; _cnt+=4){
		vr.fv = *v_vec++;
		xr.fv = *x_vec++;
		xr.fv *= vr.fv;
		yr.fv += xr.fv;
	}
	__asm__ __volatile__(" # dmv_vec end");
	return (yr.f._f0 + yr.f._f1 + yr.f._f2 + yr.f._f3);
} 
#endif

#endif /* __DMV_VEC_H__ */
