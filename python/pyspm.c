#include "Python.h"

#include <sched.h>
#include <stdlib.h>

#include "dynarray.h"
#include "matrix.h"
#include "spm_delta.h"
#include "spm_delta_cv.h"
#include "spm_mt.h"
#include "spm_delta_mt.h"
#include "spm_crs.h"
#include "spmv_loops.h"
#include "spmv_loops_mt.h"
#include "cv.h"

typedef struct {
	unsigned long nnz, ncols, nrows;
	unsigned char *ctl;
	void *values;
} spm_delta_t;

typedef struct {
	PyObject_HEAD
	spm_delta_t    delta;
	unsigned long  ctl_size;
	dynarray_t     *ctl_da, *val_da;
	unsigned int   el_type;
} pyspm_delta_t;


enum {
	PYSPM_TYPE_INVALID=0,
	PYSPM_TYPE_FLOAT,
	PYSPM_TYPE_DOUBLE
};


#define pyspm_type_unknown(type) \
	(type < PYSPM_TYPE_FLOAT) || (type > PYSPM_TYPE_DOUBLE)

static char *pyspm_eltype_fmt[] = {
	[PYSPM_TYPE_FLOAT] = "f",
	[PYSPM_TYPE_DOUBLE] = "d",
};

enum {
	PYSPM_CRS_INVALID = 0,
	PYSPM_CRS_32,
	PYSPM_CRS_64
};


#define value_error(str) PyErr_SetString(PyExc_ValueError, str)
#define index_error(str) PyErr_SetString(PyExc_IndexError, str)

#define value_error_on(cond,ret)                 \
do {                                             \
	if ( cond ) {                            \
		value_error("error on " # cond); \
		return ret;                      \
	}                                        \
} while (0)

#define index_error_on(cond,ret)                 \
do {                                             \
	if ( cond ) {                            \
		value_error("error on " # cond); \
		return ret;                      \
	}                                        \
} while (0)

#define LONG_BYTES (sizeof(unsigned long)/8)


/*
 * allocate a new pyspm_delta_object
 */
static PyObject *
pyspm_delta_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
	pyspm_delta_t *self;

	self = (pyspm_delta_t *)type->tp_alloc(type, 0);
	if ( self != NULL ){
		self->delta.nnz = self->delta.ncols = self->delta.nrows = 0;
		self->delta.values = NULL;
		self->delta.ctl = NULL;
		self->ctl_size = 0;
		self->ctl_da = self->val_da = NULL;
		self->el_type = PYSPM_TYPE_DOUBLE;
	}

	return (PyObject *)self;
}

static int
pyspm_delta_init(pyspm_delta_t *self, PyObject *args, PyObject *kwds)
{
	static char *kw[] = {"type", NULL};

	if (!PyArg_ParseTupleAndKeywords(args,kwds,"|I",kw,&self->el_type)){
		return -1;
	}

	value_error_on(pyspm_type_unknown(self->el_type), -1);

	return 0;
}
PyDoc_STRVAR(delta_init_doc,
"delta structure\n\
  initialization arguments:\n\
   type=TYPE_DOUBLE : type of elements (TYPE_DOUBLE or TYPE_FLOAT)");

static void
pyspm_delta_dealloc(PyObject *obj)
{
	pyspm_delta_t *self = (pyspm_delta_t *)obj;

	if (self->ctl_da )
		self->delta.ctl = dynarray_destroy(self->ctl_da);

	if (self->val_da)
		self->delta.values = dynarray_destroy(self->val_da);

	if (self->delta.ctl)
		free(self->delta.ctl);

	if (self->delta.values)
		free(self->delta.values);

	obj->ob_type->tp_free(obj);
}

static PyObject *
pyspm_delta_string(pyspm_delta_t *self)
{
	return PyString_FromFormat("nnz:%lu ncols:%lu nrows:%lu "
	                           "values:%p ctl:%p "
				   "val_da:%p ctl_da:%p",
				    self->delta.nnz,
				    self->delta.ncols, self->delta.nrows,
				    self->delta.values, self->delta.ctl,
				    self->ctl_da, self->val_da);
}
PyDoc_STRVAR(delta_string_doc, "string() : return some debugging info");

static PyObject *
pyspm_delta_prepare(pyspm_delta_t *self)
{
	value_error_on(self->delta.values, NULL);
	value_error_on(self->ctl_da, NULL);

	/* XXX: dynarray will call exit() if memory allocation fails */
	self->ctl_da = dynarray_create(sizeof(unsigned char), 4096);

	if ( self->el_type == PYSPM_TYPE_FLOAT ){
		self->val_da = dynarray_create(sizeof(float), 4096);
	} else if ( self->el_type == PYSPM_TYPE_DOUBLE) {
		self->val_da = dynarray_create(sizeof(double), 4096);
	}

	Py_RETURN_NONE;
}
PyDoc_STRVAR(delta_prepare_doc, "prepare a delta for parsing");

static PyObject *
pyspm_delta_addctl(pyspm_delta_t *self, PyObject *args)
{
	unsigned char ctl_val, *ctl;

	if ( !PyArg_ParseTuple(args, "B", &ctl_val) ){
		return NULL;
	}

	value_error_on(!self->ctl_da, NULL);

	ctl = dynarray_alloc(self->ctl_da);
	*ctl = ctl_val;

	Py_RETURN_NONE;
}
PyDoc_STRVAR(delta_addctl_doc, "add a control item");

static PyObject *
pyspm_delta_de_add_unit(pyspm_delta_t *self, PyObject *args, PyObject *kwds)
{
	value_error_on(!self->ctl_da, NULL);

	unsigned int size=0, nr_delta=0;
	unsigned long ci_delta=0;
	static char *kw[] = {"size", "ci_delta", "nr_delta", NULL};
	static char fmt[] = "Ik|I";

	if (!PyArg_ParseTupleAndKeywords(args, kwds, fmt, kw,
	                                 &size, &ci_delta, &nr_delta)){
		return NULL;
	}

	//printf("size=%lu max_size=%lu\n", (unsigned long)size, (unsigned long)SPM_DELTA_CTL_SIZE_MAX);
	value_error_on(size > SPM_DELTA_CTL_SIZE_MAX, NULL);
	//value_error_on(ci_delta == 0, NULL);

	unsigned char *ctl = dynarray_alloc_nr(self->ctl_da, 2);
	unsigned char *ctl_flags = &ctl[0];
	unsigned char *ctl_size = &ctl[1];

	*ctl_flags = 0;
	*ctl_size = (unsigned char)size;

	if (nr_delta > 1){
		spm_delta_fl_setnr(ctl_flags);
		da_uc_put_ul(self->ctl_da, nr_delta);
	} else if (nr_delta == 1){
		spm_delta_fl_setnrseq(ctl_flags);
		spm_delta_fl_setnr(ctl_flags);
	}

	da_uc_put_ul(self->ctl_da, ci_delta);

	//printf("%s: flags:%u (%p) size:%u (%p) \n", __FUNCTION__, *ctl_flags, ctl_flags, *ctl_size, ctl_size);

	Py_RETURN_NONE;
}
PyDoc_STRVAR(delta_de_add_unit_doc, "add control information for a dense unit");

static PyObject *
pyspm_delta_sp_add_unit_jmp(pyspm_delta_t *self, PyObject *args, PyObject *kwds)
{
	value_error_on(!self->ctl_da, NULL);

	unsigned int nr_delta=0, size=0;
	unsigned long ci_delta=0, ci_delta_max=0;
	static char *kw[] = {"size", "ci_delta", "ci_delta_max", "nr_delta", NULL};
	static char *fmt = "Ikk|I";

	if (!PyArg_ParseTupleAndKeywords(args, kwds, fmt, kw,
	                                 &size, &ci_delta, &ci_delta_max, &nr_delta)){
		return NULL;
	}

	value_error_on(size > SPM_DELTA_CTL_SIZE_MAX, NULL);
	//value_error_on(ci_delta_max == 0, NULL);

	unsigned char ci_cisize = spm_delta_cisize(ci_delta_max);

	unsigned char *ctl = dynarray_alloc_nr(self->ctl_da, 2);
	unsigned char *ctl_flags = &ctl[0];
	unsigned char *ctl_size = &ctl[1];

	*ctl_flags = 0;
	spm_delta_fl_setsp(ctl_flags);
	spm_delta_fl_setcisize(ctl_flags, ci_cisize);
	*ctl_size = (unsigned char)size;

	if (nr_delta > 1){
		spm_delta_fl_setnr(ctl_flags);
		da_uc_put_ul(self->ctl_da, nr_delta);
	} else if (nr_delta == 1){
		spm_delta_fl_setnrseq(ctl_flags);
		spm_delta_fl_setnr(ctl_flags);
	}

	da_uc_put_ul(self->ctl_da, ci_delta);

	Py_RETURN_NONE;
}

static PyObject *
pyspm_delta_sp_add_unit(pyspm_delta_t *self, PyObject *args, PyObject *kwds)
{
	value_error_on(!self->ctl_da, NULL);

	unsigned int nr_delta=0, ci_delta_max=0, size=0;
	static char *kw[] = {"size", "ci_delta_max", "nr_delta", NULL};
	static char *fmt = "Ik|I";

	if (!PyArg_ParseTupleAndKeywords(args, kwds, fmt, kw,
	                                 &size, &ci_delta_max, &nr_delta)){
		return NULL;
	}

	value_error_on(size > SPM_DELTA_CTL_SIZE_MAX, NULL);
	//value_error_on(ci_delta_max == 0, NULL);

	unsigned char ci_cisize = spm_delta_cisize(ci_delta_max);

	unsigned char *ctl = dynarray_alloc_nr(self->ctl_da, 2);
	unsigned char *ctl_flags = &ctl[0];
	unsigned char *ctl_size = &ctl[1];

	*ctl_flags = 0;
	spm_delta_fl_setsp(ctl_flags);
	spm_delta_fl_setcisize(ctl_flags, ci_cisize);
	*ctl_size = (unsigned char)size;

	if (nr_delta > 1){
		spm_delta_fl_setnr(ctl_flags);
		da_uc_put_ul(self->ctl_da, nr_delta);
	} else if (nr_delta == 1){
		spm_delta_fl_setnrseq(ctl_flags);
		spm_delta_fl_setnr(ctl_flags);
	}

	Py_RETURN_NONE;
}


static PyObject *
pyspm_delta_sp_add_cols(pyspm_delta_t *self, PyObject *args, PyObject *kwds)
{
	unsigned long size, i;
	unsigned long ci_delta_max=0;
	PyObject *ci_delta_l = NULL;

	static char *kw[] = {"ci_delta_max", "ci_delta_l", NULL};
	static char *fmt = "kO!";

	value_error_on(!self->ctl_da, NULL);

	if ( !PyArg_ParseTupleAndKeywords(args, kwds, fmt, kw,
	                                  &ci_delta_max,
	                                  &PyList_Type, &ci_delta_l)) {
		return NULL;
	}

	value_error_on(ci_delta_l == NULL, NULL);
	unsigned char ci_cisize = spm_delta_cisize(ci_delta_max);
	unsigned long ci_size = spm_delta_cisize_size(ci_cisize);

	size = PyList_Size(ci_delta_l);
	if ( !size ){
		Py_RETURN_NONE; // Support jmp semantics
	}
	//value_error_on(size == 0, NULL);

	unsigned long *ci_deltas_from = malloc(sizeof(unsigned long)*size);
	if ( !ci_deltas_from ){
		PyErr_SetString(PyExc_MemoryError, "ci_deltas allocation");
		return NULL;
	}

	for ( i=0; i<size; i++){
		PyObject *py_ci_delta = PyList_GetItem(ci_delta_l, i);
		if ( !PyInt_Check(py_ci_delta) ){
			value_error("ci_delta_l value not int");
			free(ci_deltas_from);
			return NULL;
		}

		unsigned long ci_delta = PyInt_AsUnsignedLongMask(py_ci_delta);
		if ( ci_delta > ci_delta_max ){
			value_error("ci_delta_max is not correct");
			free(ci_deltas_from);
			return NULL;
		}

		ci_deltas_from[i] = ci_delta;
	}

	void *ci_deltas_to = dynarray_alloc_nr_aligned(self->ctl_da,size*ci_size, ci_size);
	spm_delta_cisize_copy(ci_deltas_to, ci_deltas_from, size, ci_cisize);

	free(ci_deltas_from);
	Py_RETURN_NONE;
}

static PyObject *
pyspm_delta_getctl(pyspm_delta_t *self, PyObject *args)
{
	unsigned long idx;
	unsigned char ret;

	if ( !PyArg_ParseTuple(args, "k", &idx) ){
		return NULL;
	}

	if ( self->ctl_da ){
		index_error_on(idx >= dynarray_size(self->ctl_da), NULL);
		ret = *((unsigned char *)dynarray_get(self->ctl_da, idx));
	} else if ( self->delta.ctl ) {
		index_error_on(idx >= self->ctl_size, NULL);
		ret = self->delta.ctl[idx];
	} else {
		index_error_on(1, NULL);
	}

	return PyInt_FromLong((long)ret);
}

static PyObject *
pyspm_delta_addval(pyspm_delta_t *self, PyObject *args)
{
	value_error_on(!self->val_da, NULL);

	void *elem = dynarray_alloc(self->val_da);
	if ( !PyArg_ParseTuple(args, pyspm_eltype_fmt[self->el_type], elem) ){
		dynarray_dealloc(self->ctl_da);
		return NULL;
	}

	self->delta.nnz++;

	Py_RETURN_NONE;
}

static PyObject *
pyspm_delta_addvals(pyspm_delta_t *self, PyObject *args, PyObject *kwds)
{
	PyObject *values_l = NULL;
	static char *kw[] = {"values_l", NULL};
	static char *fmt = "O!";
	unsigned long size, i;

	value_error_on(!self->val_da, NULL);

	if ( !PyArg_ParseTupleAndKeywords(args, kwds, fmt, kw,
	                                  &PyList_Type, &values_l) ) {
		return NULL;
	}

	size = PyList_Size(values_l);
	value_error_on(size == 0, NULL);
	void *elems = dynarray_alloc_nr(self->val_da, size);


	for ( i=0 ; i < size ; i++){
		PyObject *py_val = PyList_GetItem(values_l, i);
		double val;
		if ( !PyFloat_Check(py_val) ) {
			value_error("values_l value not float");
			dynarray_dealloc_nr(self->val_da, size);
			return NULL;
		}

		val = PyFloat_AsDouble(py_val);
		//printf("inserting val=%lf\n", val);
		if ( self->el_type == PYSPM_TYPE_FLOAT){
			float *_val = (float *)elems;
			*_val++ = (float)val;
			elems = _val;
		} else {
			double *_val = (double *)elems;
			*_val++ = val;
			elems = _val;
		}
	}

	self->delta.nnz += size;
	Py_RETURN_NONE;
}

static PyObject *
pyspm_delta_getval(pyspm_delta_t *self, PyObject *args)
{
	unsigned long idx;
	double ret;

	if ( !PyArg_ParseTuple(args, "k", &idx) ){
		return NULL;
	}

	if ( self->val_da ){
		index_error_on( idx >= dynarray_size(self->val_da), NULL);
		if ( self->el_type == PYSPM_TYPE_FLOAT){
			float _ret = *((float *)dynarray_get(self->val_da, idx));
			ret = (double)_ret;
		} else {
			ret = *(double *)dynarray_get(self->val_da, idx);
		}
	} else if ( self->delta.values ) {
		index_error_on( idx >= self->delta.nnz, NULL);
		if ( self->el_type == PYSPM_TYPE_FLOAT){
			float _ret = *((float *)self->delta.values + idx);
			ret = (double)_ret;
		} else {
			ret = *((double *)self->delta.values + idx);
		}
	} else {
		index_error_on(1, NULL);
	}

	return PyFloat_FromDouble(ret);
}

#if 0
static PyObject *pyspm_delta_popctl(pyspm_delta_t *self)
{
	Py_RETURN_NONE;
}

static PyObject *pyspm_delta_popval(pyspm_delta_t *self)
{
	Py_RETURN_NONE;
}
#endif

static PyObject *
pyspm_delta_finalize(pyspm_delta_t *self, PyObject *args)
{
	unsigned long nnz, nrows, ncols, ctl_size;

	if ( !PyArg_ParseTuple(args, "kkk", &nnz, &nrows, &ncols) ){
		return NULL;
	}

	ctl_size = dynarray_size(self->ctl_da);
	self->delta.ctl = dynarray_destroy(self->ctl_da);
	self->delta.values = dynarray_destroy(self->val_da);
	self->ctl_da = self->val_da = NULL;

	//printf("nnz=%lu, self->delta.nnz=%lu ctl_size:%lu\n", nnz, self->delta.nnz, ctl_size);
	value_error_on(nnz != self->delta.nnz, NULL);

	self->delta.nnz = nnz;
	self->delta.nrows = nrows;
	self->delta.ncols = ncols;
	self->ctl_size = ctl_size;

	Py_RETURN_NONE;
}

static PyObject *
pyspm_delta_getprops(pyspm_delta_t *self)
{
	PyObject *d, *item;

	if ( !(d = PyDict_New()) ){
		return NULL;
	}

	item = PyLong_FromUnsignedLong(self->delta.nnz);
	PyDict_SetItemString(d, "nnz", item);
	Py_DECREF(item);

	item = PyLong_FromUnsignedLong(self->delta.ncols);
	PyDict_SetItemString(d, "ncols", item);
	Py_DECREF(item);

	item = PyLong_FromUnsignedLong(self->delta.nrows);
	PyDict_SetItemString(d, "nrows", item);
	Py_DECREF(item);

	item = PyLong_FromUnsignedLong(self->ctl_size);
	PyDict_SetItemString(d, "ctl_size", item);
	Py_DECREF(item);

	return d;
}

static PyObject *
pyspm_delta_check(pyspm_delta_t *self, PyObject *args)
{
	char *file=NULL;
	unsigned long loops=1;

	if ( !PyArg_ParseTuple(args, "s|k", &file, &loops) ){
		return NULL;
	}

	value_error_on(file == NULL, NULL);

	unsigned long rows_nr, cols_nr, nz_nr;
	if ( self->el_type == PYSPM_TYPE_FLOAT){
		spm_crs64_float_t *crs;
		crs = spm_crs64_float_init_mmf(file, &rows_nr, &cols_nr, &nz_nr);
		spmv_float_check_loop(
			crs, &self->delta,
			spm_crs64_float_multiply, spm_delta_float_multiply,
			loops, cols_nr, nz_nr);
		spm_crs64_float_destroy(crs);
	} else {
		spm_crs64_double_t *crs;
		crs = spm_crs64_double_init_mmf(file, &rows_nr, &cols_nr, &nz_nr);
		spmv_double_check_loop(
			crs, &self->delta,
			spm_crs64_double_multiply, spm_delta_double_multiply,
			loops, cols_nr, nz_nr);
		spm_crs64_double_destroy(crs);
	}

	Py_RETURN_NONE;
}

static PyObject *
pyspm_delta_jmp_check(pyspm_delta_t *self, PyObject *args)
{
	char *file=NULL;
	unsigned long loops=1;

	if ( !PyArg_ParseTuple(args, "s|k", &file, &loops) ){
		return NULL;
	}

	value_error_on(file == NULL, NULL);

	unsigned long rows_nr, cols_nr, nz_nr;
	if ( self->el_type == PYSPM_TYPE_FLOAT){
		spm_crs64_float_t *crs;
		crs = spm_crs64_float_init_mmf(file, &rows_nr, &cols_nr, &nz_nr);
		spmv_float_check_loop(
			crs, &self->delta,
			spm_crs64_float_multiply, spm_delta_float_jmp_multiply,
			loops, cols_nr, nz_nr);
		spm_crs64_float_destroy(crs);
	} else {
		spm_crs64_double_t *crs;
		crs = spm_crs64_double_init_mmf(file, &rows_nr, &cols_nr, &nz_nr);
		spmv_double_check_loop(
			crs, &self->delta,
			spm_crs64_double_multiply, spm_delta_double_jmp_multiply,
			loops, cols_nr, nz_nr);
		spm_crs64_double_destroy(crs);
	}

	Py_RETURN_NONE;
}

static PyObject *
pyspm_delta_bench(pyspm_delta_t *self, PyObject *args)
{
	unsigned long loops=128;
	double ret;

	if ( !PyArg_ParseTuple(args, "|k", &loops) ){
		return NULL;
	}

	if ( self->el_type == PYSPM_TYPE_FLOAT){
		ret = spmv_float_bench_loop(spm_delta_float_multiply,
		                            &self->delta, loops, self->delta.ncols);
	} else {
		ret = spmv_double_bench_loop(spm_delta_double_multiply,
		                             &self->delta, loops, self->delta.ncols);
	}

	return Py_BuildValue("(dd)", ret, (loops*self->delta.nnz*2)/(1000*1000*ret));
}

static PyObject *
pyspm_delta_jmp_bench(pyspm_delta_t *self, PyObject *args)
{
	unsigned long loops=128;
	double ret;

	if ( !PyArg_ParseTuple(args, "|k", &loops) ){
		return NULL;
	}

	if ( self->el_type == PYSPM_TYPE_FLOAT){
		ret = spmv_float_bench_loop(spm_delta_float_jmp_multiply,
		                            &self->delta, loops, self->delta.ncols);
	} else {
		ret = spmv_double_bench_loop(spm_delta_double_jmp_multiply,
		                             &self->delta, loops, self->delta.ncols);
	}

	return Py_BuildValue("(dd)", ret, (loops*self->delta.nnz*2)/(1000*1000*ret));
}

static PyObject *
pyspm_delta_jmp_mt_bench(pyspm_delta_t *self, PyObject *args)
{
	unsigned long loops=128;
	double ret;


	if ( !PyArg_ParseTuple(args, "|k", &loops) ){
		return NULL;
	}

	spm_mt_t *spm_mt = spm_delta_jmp_mt_partition(self->delta.ctl,
	                                          self->delta.nnz,
	                                          &self->delta);

	if ( self->el_type == PYSPM_TYPE_FLOAT){
		ret = spmv_float_bench_mt_loop(spm_delta_mt_float_jmp_multiply,
		                             spm_mt,
		                             loops,
					     self->delta.ncols);
	} else {
		ret = spmv_double_bench_mt_loop(spm_delta_mt_double_jmp_multiply,
		                             spm_mt,
		                             loops,
					     self->delta.ncols);
	}

	return Py_BuildValue("(dd)", ret, (loops*self->delta.nnz*2)/(1000*1000*ret));
}

static PyObject *
pyspm_delta_jmp_mt_check(pyspm_delta_t *self, PyObject *args)
{
	char *file=NULL;
	unsigned long loops=1;


	if ( !PyArg_ParseTuple(args, "s|k", &file, &loops) ){
		return NULL;
	}

	value_error_on(file==NULL, NULL);

	spm_mt_t *spm_mt = spm_delta_jmp_mt_partition(self->delta.ctl,
	                                          self->delta.nnz,
	                                          &self->delta);

	unsigned long rows_nr, cols_nr, nz_nr;
	if ( self->el_type == PYSPM_TYPE_FLOAT){
		spm_crs32_float_t *crs;
		crs = spm_crs32_float_init_mmf(file, &rows_nr, &cols_nr, &nz_nr);
		spmv_float_check_mt_loop(crs, spm_mt,
		                         spm_crs32_float_multiply,
		                         spm_delta_mt_float_jmp_multiply,
		                         loops, self->delta.ncols);
	} else {
		spm_crs32_double_t *crs;
		crs = spm_crs32_double_init_mmf(file, &rows_nr, &cols_nr, &nz_nr);
		spmv_double_check_mt_loop(crs, spm_mt,
		                          spm_crs32_double_multiply,
		                          spm_delta_mt_double_jmp_multiply,
		                          loops, self->delta.ncols);
	}

	Py_RETURN_NONE;
}

/*
 * pyspm.delta methods
 */
static PyMethodDef pyspm_delta_methods[] = {
	{"string", (PyCFunction)pyspm_delta_string, METH_NOARGS, delta_string_doc},
	{"getprops", (PyCFunction)pyspm_delta_getprops, METH_NOARGS, ""},
	{"prepare", (PyCFunction)pyspm_delta_prepare, METH_NOARGS, delta_prepare_doc},
	{"addctl", (PyCFunction)pyspm_delta_addctl, METH_VARARGS, delta_addctl_doc },
	{"addval", (PyCFunction)pyspm_delta_addval, METH_VARARGS, "" },
	{"addvals", (PyCFunction)pyspm_delta_addvals, METH_KEYWORDS, "" },
	{"add_de_unit", (PyCFunction)pyspm_delta_de_add_unit, METH_KEYWORDS, delta_de_add_unit_doc},
	{"add_sp_unit", (PyCFunction)pyspm_delta_sp_add_unit, METH_KEYWORDS, "" },
	{"add_sp_unit_jmp", (PyCFunction)pyspm_delta_sp_add_unit_jmp, METH_KEYWORDS, "" },
	{"add_sp_cols", (PyCFunction)pyspm_delta_sp_add_cols, METH_KEYWORDS, "" },
	{"finalize", (PyCFunction)pyspm_delta_finalize, METH_VARARGS, "" },
	{"getctl", (PyCFunction)pyspm_delta_getctl, METH_VARARGS, "" },
	{"getval", (PyCFunction)pyspm_delta_getval, METH_VARARGS, "" },
	{"check", (PyCFunction)pyspm_delta_check, METH_VARARGS, "" },
	{"check_jmp", (PyCFunction)pyspm_delta_jmp_check, METH_VARARGS, "" },
	{"bench", (PyCFunction)pyspm_delta_bench, METH_VARARGS, "" },
	{"bench_jmp", (PyCFunction)pyspm_delta_jmp_bench, METH_VARARGS, "" },
	{"mt_bench_jmp", (PyCFunction)pyspm_delta_jmp_mt_bench, METH_VARARGS, "" },
	{"mt_check_jmp", (PyCFunction)pyspm_delta_jmp_mt_check, METH_VARARGS, "" },
	{NULL}
};


/*
 * pyspm.delta type
 */
static PyTypeObject pyspm_delta_type = {
	PyObject_HEAD_INIT(NULL)
	0,                         /*ob_size*/
	"pyspm.delta",             /*tp_name*/
	sizeof(pyspm_delta_t),     /*tp_basicsize*/
	0,                         /*tp_itemsize*/
	pyspm_delta_dealloc,       /*tp_dealloc*/
	0,                         /*tp_print*/
	0,                         /*tp_getattr*/
	0,                         /*tp_setattr*/
	0,                         /*tp_compare*/
	0,                         /*tp_repr*/
	0,                         /*tp_as_number*/
	0,                         /*tp_as_sequence*/
	0,                         /*tp_as_mapping*/
	0,                         /*tp_hash */
	0,                         /*tp_call*/
	0,                         /*tp_str*/
	0,                         /*tp_getattro*/
	0,                         /*tp_setattro*/
	0,                         /*tp_as_buffer*/
	Py_TPFLAGS_DEFAULT,        /*tp_flags*/
	delta_init_doc,            /* tp_doc */
	0,                         /* tp_traverse */
	0,                         /* tp_clear */
	0,                         /* tp_richcompare */
	0,                         /* tp_weaklistoffset */
	0,                         /* tp_iter */
	0,                         /* tp_iternext */
	pyspm_delta_methods,       /* tp_methods */
	0,                         /* tp_members */
	0,                         /* tp_getset */
	0,                         /* tp_base */
	0,                         /* tp_dict */
	0,                         /* tp_descr_get */
	0,                         /* tp_descr_set */
	0,                         /* tp_dictoffset */
	(initproc)pyspm_delta_init,/* tp_init */
	0,                         /* tp_alloc */
	pyspm_delta_new,           /* tp_new */
};

typedef struct {
	unsigned long nnz, ncols, nrows;
	unsigned long ctl_size, ci_size, vi_size, v_size;
	unsigned char *ctl;
	unsigned char *c_indices;
	unsigned char *v_indices;
	void          *values;
} spm_delta_cv_t;

typedef struct {
	PyObject_HEAD
	spm_delta_cv_t delta_cv;
	dynarray_t     *ctl_da, *ci_da, *vi_da, *val_da;
	unsigned int   el_type;
} pyspm_delta_cv_t;

#define FREE_IF_NOTNULL(ptr)       \
	do {                       \
		if (ptr != NULL){  \
			free(ptr); \
		}                  \
	} while(0)

#define DA_DESTROY_IF_NOTNULL(da, ptr)              \
	do {                                        \
		if (da != NULL){                    \
			ptr = dynarray_destroy(da); \
		}                                   \
	} while(0)

static PyObject *
pyspm_delta_cv_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
	pyspm_delta_cv_t *self;

	self = (pyspm_delta_cv_t *)type->tp_alloc(type, 0);
	if ( self != NULL ){
		self->el_type = PYSPM_TYPE_DOUBLE;
		spm_delta_cv_t *d = &self->delta_cv;
		d->nnz = d->ncols = d->nrows = 0;
		d->values = d->c_indices = d->v_indices = d->ctl = NULL;
	}

	return (PyObject *)self;
}

static PyObject *
pyspm_delta_cv_init(pyspm_delta_cv_t *self, PyObject *args, PyObject *kwds)
{
	return 0;
}

static void
pyspm_delta_cv_dealloc(PyObject *obj)
{
	pyspm_delta_cv_t *self = (pyspm_delta_cv_t *)obj;

	DA_DESTROY_IF_NOTNULL(self->ctl_da, self->delta_cv.ctl);
	DA_DESTROY_IF_NOTNULL(self->ci_da, self->delta_cv.c_indices);
	DA_DESTROY_IF_NOTNULL(self->vi_da, self->delta_cv.v_indices);
	DA_DESTROY_IF_NOTNULL(self->val_da, self->delta_cv.values);

	FREE_IF_NOTNULL(self->delta_cv.ctl);
	FREE_IF_NOTNULL(self->delta_cv.c_indices);
	FREE_IF_NOTNULL(self->delta_cv.v_indices);
	FREE_IF_NOTNULL(self->delta_cv.values);

	obj->ob_type->tp_free(obj);
}

static PyObject *
pyspm_delta_cv_prepare(pyspm_delta_cv_t *self)
{
	value_error_on(self->ctl_da, NULL);
	value_error_on(self->ci_da, NULL);
	value_error_on(self->vi_da, NULL);
	value_error_on(self->val_da, NULL);

	value_error_on(self->delta_cv.ctl, NULL);
	value_error_on(self->delta_cv.c_indices, NULL);
	value_error_on(self->delta_cv.v_indices, NULL);
	value_error_on(self->delta_cv.values, NULL);
	value_error_on(pyspm_type_unknown(self->el_type), NULL);

	self->ctl_da = dynarray_create(sizeof(unsigned char), 4096);
	self->ci_da = dynarray_create(sizeof(unsigned char), 4096);
	self->vi_da = dynarray_create(sizeof(unsigned char), 4096);


	if ( self->el_type == PYSPM_TYPE_FLOAT ){
		self->val_da = dynarray_create(sizeof(float), 4096);
	} else if ( self->el_type == PYSPM_TYPE_DOUBLE) {
		self->val_da = dynarray_create(sizeof(double), 4096);
	}

	Py_RETURN_NONE;
}

static PyObject *
pyspm_delta_cv_add_unit(pyspm_delta_cv_t *self, PyObject *args, PyObject *kwds)
{
	unsigned int u_s=0, ci_s=0, vi_s=0, sp=1;
	unsigned long ci_jmp=0, ri_jmp=0;

	static char *kw[] = {"u_s","ci_s","vi_s", "ci_jmp", "ri_jmp","sp",NULL};
	static char *fmt = "IIIkkI";

	if ( !PyArg_ParseTupleAndKeywords(args, kwds, fmt, kw,
	                                  &u_s, &ci_s, &vi_s,
	                                  &ci_jmp, &ri_jmp,
                                          &sp) ){
		return NULL;
	}

	value_error_on(u_s > SPM_DELTA_CTL_SIZE_MAX, NULL);

	unsigned char *ctl = dynarray_alloc_nr(self->ctl_da, 2);
	unsigned char *ctl_flags = &ctl[0];
	unsigned char *ctl_size = &ctl[1];

	*ctl_flags = 0;
	*ctl_size = (unsigned char)u_s;
	if ( sp == 1){
		spm_delta_cv_fl_set_cvi_sp(ctl_flags, ci_s, vi_s);
	} else {
		spm_delta_cv_fl_set_cvi_de(ctl_flags, vi_s);
	}

	if (ri_jmp > 1){
		spm_delta_fl_setnr(ctl_flags);
		da_uc_put_ul(self->ctl_da, ri_jmp);
	} else if (ri_jmp == 1){
		spm_delta_fl_setnrseq(ctl_flags);
		spm_delta_fl_setnr(ctl_flags);
	}

	da_uc_put_ul(self->ctl_da, ci_jmp);

	Py_RETURN_NONE;
}

static int
pyspm_delta_add_uindices(dynarray_t *da, PyObject *i_list, int i_size)
{
	int i;
	void *index_to;
	unsigned int i_bytes = spm_delta_cisize_size(i_size);
	unsigned long size = PyList_Size(i_list);
	if (size == 0){
		return 0;
	}

	index_to = dynarray_alloc_nr(da, size*i_bytes);

	for (i=0; i<size; i++){
		PyObject *py_i = PyList_GetItem(i_list, i);
		if ( !PyInt_Check(py_i) ){
			value_error("not int");
			return -1;
		}

		unsigned long val_i = PyInt_AsUnsignedLongMask(py_i);
		spm_delta_cisize_set(index_to, val_i, i_size);
		index_to += i_bytes;
	}

	return 0;
}

static PyObject *
pyspm_delta_cv_add_ci(pyspm_delta_cv_t *self, PyObject *args, PyObject *kwds)
{
	int ci_size=-1, err;
	PyObject  *ci_delta_l = NULL;

	static char *kw[] = {"ci_s", "ci_delta_l", NULL};
	static char *fmt = "iO!";

	value_error_on(!self->ci_da, NULL);

	if ( !PyArg_ParseTupleAndKeywords(args, kwds, fmt, kw,
	                                  &ci_size,
					  &PyList_Type, &ci_delta_l)) {
		return NULL;
	}
	value_error_on(ci_delta_l == NULL, NULL);

	err = pyspm_delta_add_uindices(self->ci_da, ci_delta_l, ci_size);
	if ( err ){
		return NULL;
	}

	Py_RETURN_NONE;
}

#if 0
static int
pyspm_delta_add_sindices(dynarray_t *da, PyObject *i_list, int i_size)
{
	int i;
	void *index_to;
	unsigned int i_bytes = spm_delta_cisize_size(i_size);
	unsigned long size = PyList_Size(i_list);
	if (size == 0){
		return 0;
	}

	index_to = dynarray_alloc_nr(da, size*i_bytes);

	for (i=0; i<size; i++){
		PyObject *py_i = PyList_GetItem(i_list, i);
		if ( !PyInt_Check(py_i) ){
			value_error("not int");
			return -1;
		}

		long val_i = PyInt_AS_LONG(py_i);
		si_set(index_to, val_i, i_size);
		index_to += i_bytes;
	}

	return 0;
}
#endif

static PyObject *
pyspm_delta_cv_add_vi(pyspm_delta_cv_t *self, PyObject *args, PyObject *kwds)
{
	int vi_size=-1, err;
	PyObject  *vi_l = NULL;

	static char *kw[] = {"vi_s", "vi_l", NULL};
	static char *fmt = "iO!";

	value_error_on(!self->vi_da, NULL);

	if ( !PyArg_ParseTupleAndKeywords(args, kwds, fmt, kw,
	                                  &vi_size,
					  &PyList_Type, &vi_l)) {
		return NULL;
	}
	value_error_on(vi_l == NULL, NULL);
	self->delta_cv.nnz += PyList_Size(vi_l);

	err = pyspm_delta_add_uindices(self->vi_da, vi_l, vi_size);
	if ( err ){
		return NULL;
	}

	Py_RETURN_NONE;
}

static PyObject *
pyspm_delta_cv_add_vals(pyspm_delta_cv_t *self, PyObject *args, PyObject *kwds)
{
	PyObject *values_l = NULL;
	static char *kw[] = {"v_l", NULL};
	static char *fmt = "O!";

	value_error_on(!self->val_da, NULL);

	if ( !PyArg_ParseTupleAndKeywords(args, kwds, fmt, kw,
	                                  &PyList_Type, &values_l) ){
		return NULL;
	}

	unsigned long size = PyList_Size(values_l);
	value_error_on(size == 0, NULL);
	void *elems = dynarray_alloc_nr(self->val_da, size);

	int i;
	for (i=0 ; i<size ; i++){
		PyObject *py_val = PyList_GetItem(values_l, i);
		double val;
		if ( !PyFloat_Check(py_val) ) {
			value_error("values_l value not float");
			dynarray_dealloc_nr(self->val_da, size);
			return NULL;
		}

		val = PyFloat_AsDouble(py_val);
		//printf("inserting val=%lf\n", val);
		if ( self->el_type == PYSPM_TYPE_FLOAT){
			float *_val = (float *)elems;
			*_val++ = (float)val;
			elems = _val;
		} else {
			double *_val = (double *)elems;
			*_val++ = val;
			elems = _val;
		}
	}

	Py_RETURN_NONE;
}

static PyObject *
pyspm_delta_cv_finalize(pyspm_delta_cv_t *self, PyObject *args)
{
	unsigned long nnz, nrows, ncols;

	if ( !PyArg_ParseTuple(args, "kkk", &nnz, &nrows, &ncols) ){
		return NULL;
	}
	//
	//printf("nnz=%lu, self->delta.nnz=%lu\n", nnz, self->delta.nnz);
	value_error_on(nnz != self->delta_cv.nnz, NULL);

	spm_delta_cv_t *d = &self->delta_cv;

	d->ctl_size = dynarray_size(self->ctl_da);
	d->ctl = dynarray_destroy(self->ctl_da);

	d->v_size = dynarray_size(self->val_da);
	d->values = dynarray_destroy(self->val_da);

	d->ci_size = dynarray_size(self->ci_da);
	d->c_indices = dynarray_destroy(self->ci_da);

	d->vi_size = dynarray_size(self->vi_da);
	d->v_indices = dynarray_destroy(self->vi_da);

	self->vi_da = self->ci_da = self->val_da = self->ctl_da = NULL;

	self->delta_cv.nnz = nnz;
	self->delta_cv.nrows = nrows;
	self->delta_cv.ncols = ncols;

	Py_RETURN_NONE;
}

static PyObject *
pyspm_delta_cv_check(pyspm_delta_cv_t *self, PyObject *args)
{
	char *file=NULL;
	unsigned long loops=1;

	if ( !PyArg_ParseTuple(args, "s|k", &file, &loops) ){
		return NULL;
	}

	value_error_on(file == NULL, NULL);

	unsigned long rows_nr, cols_nr, nz_nr;
	if ( self->el_type == PYSPM_TYPE_FLOAT){
		spm_crs64_float_t *crs;
		crs = spm_crs64_float_init_mmf(file, &rows_nr, &cols_nr, &nz_nr);
		spmv_float_check_loop(
			crs, &self->delta_cv,
			spm_crs64_float_multiply, spm_delta_cv_float_jmp_multiply,
			loops, cols_nr, nz_nr);
		spm_crs64_float_destroy(crs);
	} else {
		spm_crs64_double_t *crs;
		crs = spm_crs64_double_init_mmf(file, &rows_nr, &cols_nr, &nz_nr);
		spmv_double_check_loop(
			crs, &self->delta_cv,
			spm_crs64_double_multiply, spm_delta_cv_double_jmp_multiply,
			loops, cols_nr, nz_nr);
		spm_crs64_double_destroy(crs);
	}

	Py_RETURN_NONE;
}

static PyObject *
pyspm_delta_cv_bench(pyspm_delta_cv_t *self, PyObject *args)
{
	unsigned long loops=128;
	double ret;

	if ( !PyArg_ParseTuple(args, "|k", &loops) ){
		return NULL;
	}

	if ( self->el_type == PYSPM_TYPE_FLOAT){
		ret = spmv_float_bench_loop(spm_delta_cv_float_jmp_multiply,
		                            &self->delta_cv, loops, self->delta_cv.ncols);
	} else {
		ret = spmv_double_bench_loop(spm_delta_cv_double_jmp_multiply,
		                             &self->delta_cv, loops, self->delta_cv.ncols);
	}

	ret = (double)(loops*self->delta_cv.nnz*2)/(1000*1000*ret);
	return PyFloat_FromDouble(ret);
}

static PyMethodDef pyspm_delta_cv_methods[] = {
	{"prepare", (PyCFunction) pyspm_delta_cv_prepare, METH_NOARGS,  ""},
	{"add_unit", (PyCFunction) pyspm_delta_cv_add_unit, METH_KEYWORDS, ""},
	{"add_ci", (PyCFunction) pyspm_delta_cv_add_ci, METH_KEYWORDS, ""},
	{"add_vi", (PyCFunction) pyspm_delta_cv_add_vi, METH_KEYWORDS, ""},
	{"add_vals", (PyCFunction) pyspm_delta_cv_add_vals, METH_KEYWORDS, ""},
	{"finalize", (PyCFunction) pyspm_delta_cv_finalize, METH_KEYWORDS, ""},
	{"check",   (PyCFunction) pyspm_delta_cv_check,   METH_VARARGS, "" },
	{"bench",   (PyCFunction) pyspm_delta_cv_bench,   METH_VARARGS, "" },
	{NULL}
};

/*
 * pyspm.delta_cv type
 */
static PyTypeObject pyspm_delta_cv_type = {
	PyObject_HEAD_INIT(NULL)
	0,                              /*ob_size*/
	"pyspm.delta_cv",               /*tp_name*/
	sizeof(pyspm_delta_cv_t),       /*tp_basicsize*/
	0,                              /*tp_itemsize*/
	pyspm_delta_cv_dealloc,         /*tp_dealloc*/
	0,                              /*tp_print*/
	0,                              /*tp_getattr*/
	0,                              /*tp_setattr*/
	0,                              /*tp_compare*/
	0,                              /*tp_repr*/
	0,                              /*tp_as_number*/
	0,                              /*tp_as_sequence*/
	0,                              /*tp_as_mapping*/
	0,                              /*tp_hash */
	0,                              /*tp_call*/
	0,                              /*tp_str*/
	0,                              /*tp_getattro*/
	0,                              /*tp_setattro*/
	0,                              /*tp_as_buffer*/
	Py_TPFLAGS_DEFAULT,             /*tp_flags*/
	"delta_cv",                     /* tp_doc */
	0,                              /* tp_traverse */
	0,                              /* tp_clear */
	0,                              /* tp_richcompare */
	0,                              /* tp_weaklistoffset */
	0,                              /* tp_iter */
	0,                              /* tp_iternext */
	pyspm_delta_cv_methods,         /* tp_methods */
	0,                              /* tp_members */
	0,                              /* tp_getset */
	0,                              /* tp_base */
	0,                              /* tp_dict */
	0,                              /* tp_descr_get */
	0,                              /* tp_descr_set */
	0,                              /* tp_dictoffset */
	(initproc)pyspm_delta_cv_init,  /* tp_init */
	0,                              /* tp_alloc */
	pyspm_delta_cv_new,             /* tp_new */
};

static PyObject *
pyspm_delta_cisize(PyObject *self, PyObject *args)
{
	unsigned long cidelta = 0;
	if ( !PyArg_ParseTuple(args, "k", &cidelta) ){
		return NULL;
	}
	//value_error_on(!cidelta, NULL);

	long ret = spm_delta_cisize(cidelta);

	return PyInt_FromLong(ret);
}

static PyObject *
pyspm_delta_visize(PyObject *self, PyObject *args)
{
	unsigned long videlta = 0;
	if ( !PyArg_ParseTuple(args, "k", &videlta) ){
		return NULL;
	}
	//value_error_on(!cidelta, NULL);

	long ret = spm_delta_visize(videlta);

	return PyInt_FromLong(ret);
}

static PyObject *
pyspm_crs_bench(PyObject *self, PyObject *args)
{
	unsigned int type=PYSPM_TYPE_INVALID, crs=PYSPM_CRS_INVALID;
	unsigned long loops=128;
	char *f=NULL;
	double ret = (double)0;

	if ( !PyArg_ParseTuple(args, "sII|k", &f, &type, &crs, &loops) ){
		return NULL;
	}

	unsigned long rows_nr, cols_nr, nz_nr;
	if (type == PYSPM_TYPE_FLOAT && crs == PYSPM_CRS_32) {
		spm_crs32_float_t *crs;
		crs = spm_crs32_float_init_mmf(f, &rows_nr, &cols_nr, &nz_nr);
		ret = spmv_float_bench_loop(spm_crs32_float_multiply, crs, loops, cols_nr);
		spm_crs32_float_destroy(crs);
	} else if (type == PYSPM_TYPE_FLOAT && crs == PYSPM_CRS_64) {
		spm_crs64_float_t *crs;
		crs = spm_crs64_float_init_mmf(f, &rows_nr, &cols_nr, &nz_nr);
		ret = spmv_float_bench_loop(spm_crs64_float_multiply, crs, loops, cols_nr);
		spm_crs64_float_destroy(crs);
	} else if (type == PYSPM_TYPE_DOUBLE && crs == PYSPM_CRS_32) {
		spm_crs32_double_t *crs;
		crs = spm_crs32_double_init_mmf(f, &rows_nr, &cols_nr, &nz_nr);
		ret = spmv_double_bench_loop(spm_crs32_double_multiply, crs, loops, cols_nr);
		spm_crs32_double_destroy(crs);
	} else if (type == PYSPM_TYPE_DOUBLE && crs == PYSPM_CRS_64) {
		spm_crs64_double_t *crs;
		crs = spm_crs64_double_init_mmf(f, &rows_nr, &cols_nr, &nz_nr);
		ret = spmv_double_bench_loop(spm_crs64_double_multiply, crs, loops, cols_nr);
		spm_crs64_double_destroy(crs);
	} else {
		value_error_on(1, NULL);
	}

	ret = (double)(loops*nz_nr*2)/(1000*1000*ret);
	return PyFloat_FromDouble(ret);
}

static PyObject *
pyspm_setaffinity(PyObject *self, PyObject *args)
{
	int err;
	unsigned int cpu;
	if (!PyArg_ParseTuple(args, "I", &cpu)){
		return NULL;
	}

	cpu_set_t cpu_mask;
	CPU_ZERO(&cpu_mask);
	CPU_SET(cpu, &cpu_mask);

	err = sched_setaffinity(0, sizeof(cpu_set_t), &cpu_mask);
	if ( err ){
		perror("sched_setaffinity");
		exit(1);
	}

	Py_RETURN_NONE;
}

static PyMethodDef pyspm_methods[] = {
	{"crs_bench", pyspm_crs_bench, METH_VARARGS, "" },
	{"setaffinity", pyspm_setaffinity, METH_VARARGS, ""},
	{"delta_cisize", pyspm_delta_cisize, METH_VARARGS, ""},
	{"delta_visize", pyspm_delta_visize, METH_VARARGS, ""},
	{NULL, NULL}
};

/*
 * initialization function for pyspm module
 */
PyMODINIT_FUNC
initpyspm(void)
{
	PyObject *m;

	if (PyType_Ready(&pyspm_delta_type) < 0){
		return;
	}

	if (PyType_Ready(&pyspm_delta_cv_type) < 0){
		return;
	}

	m = Py_InitModule3("pyspm", pyspm_methods, "pyspm");

	Py_INCREF(&pyspm_delta_type);
	PyModule_AddObject(m, "delta", (PyObject *)&pyspm_delta_type);

	Py_INCREF(&pyspm_delta_cv_type);
	PyModule_AddObject(m, "delta_cv", (PyObject *)&pyspm_delta_cv_type);

	PyModule_AddIntConstant(m,"TYPE_FLOAT",PYSPM_TYPE_FLOAT);
	PyModule_AddIntConstant(m,"TYPE_DOUBLE",PYSPM_TYPE_DOUBLE);
	PyModule_AddIntConstant(m,"CRS32" ,PYSPM_CRS_32);
	PyModule_AddIntConstant(m,"CRS64" ,PYSPM_CRS_64);

	PyModule_AddIntConstant(m,"DELTA_CTL_SIZE_MAX",SPM_DELTA_CTL_SIZE_MAX);

	PyModule_AddIntConstant(m,"DELTA_CISIZE_U8", SPM_DELTA_CISIZE_U8);
	PyModule_AddIntConstant(m,"DELTA_CISIZE_U16", SPM_DELTA_CISIZE_U16);
	PyModule_AddIntConstant(m,"DELTA_CISIZE_U32", SPM_DELTA_CISIZE_U32);
	PyModule_AddIntConstant(m,"DELTA_CISIZE_U64", SPM_DELTA_CISIZE_U64);
}
