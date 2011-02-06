#ifndef __ELEM_H__
#define __ELEM_H__

enum {
	ELEM_FLOAT = 1,
	ELEM_DOUBLE
};

static inline int elem_size(unsigned type)
{
	switch (type){
		case ELEM_FLOAT:
		return sizeof(float);

		case ELEM_DOUBLE:
		return sizeof(double);

		default:
		fprintf(stderr, "%s: uknown type: %u\n", __FUNCTION__, type);
		exit(1);
	}
}

static inline void elem_set(void *dst, unsigned long dst_idx,
                            void *src, unsigned long src_idx, 
			    unsigned type)
{
	switch (type){
		case ELEM_FLOAT:
		{
			float *_src = (float *)src;
			float *_dst = (float *)dst;
			_dst[dst_idx] = _src[src_idx];
		}
		break;

		case ELEM_DOUBLE:
		{
			double *_src = (double *)src;
			double *_dst = (double *)dst;
			_dst[dst_idx] = _src[src_idx];
		}

		default:
		fprintf(stderr, "%s: uknown type: %u\n", __FUNCTION__, type);
		exit(1);
	}
}

static inline void elem_set_from_double(void *dst, unsigned long dst_idx,
                                        double val, unsigned type)
{
	switch (type){
		case ELEM_FLOAT:
		{
			float *_dst = (float *)dst;
			_dst[dst_idx] = (float)val;
		}
		break;

		case ELEM_DOUBLE:
		{
			double *_dst = (double *)dst;
			_dst[dst_idx] = val;
		}

		default:
		fprintf(stderr, "%s: uknown type: %u\n", __FUNCTION__, type);
		exit(1);
	}
}

#endif
