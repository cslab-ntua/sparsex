#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <inttypes.h>
#include <assert.h>

#include "dynarray.h"
#include "mmf.h"
#include "drle.h"

#define STATS_VAL_NR 256
#define FREQ_LIMIT 4

int point_cmp_fn(const void *o0, const void *o1)
{
	const point_t *p0, *p1;
	int64_t ret;

	p0 = o0;
	p1 = o1;

	ret = p0->y - p1->y;
	if (ret == 0){
		ret = p0->x - p1->x;
	}

	if (ret > 0){
		return 1;
	} else if (ret < 0){
		return -1;
	} else {
		return 0;
	}
}

#if 0
static void print_points(point_t *points, uint64_t points_nr)
{
	uint64_t i;
	for (i=0; i<points_nr; i++){
		printf("(%lu %lu) ", points[i].y, points[i].x);
	}
	printf("\n");
}
#endif

struct rle_stat {
	unsigned long count;
	unsigned long total;
};
typedef struct rle_stat rle_stat_t;

static void delta_encode(uint64_t *input, uint64_t *deltas, uint64_t size)
{
	uint64_t i;
	uint64_t prev = input[0];
	deltas[0] = prev;
	for (i=1; i<size; i++){
		uint64_t curr = input[i];
		deltas[i] = curr - prev;
		prev = curr;
	}
	//printf("prev=%lu input=%lu\n", deltas[0], input[0]);
}

static void print_stats(rle_stat_t *stats, unsigned long nnz)
{
	int i;
	rle_stat_t *s;
	double dnnz = (double)nnz;

	for (i=0; i<STATS_VAL_NR; i++){
		s = &stats[i];
		if (s->count == 0){
			continue;
		}
		double p = (100.0*(double)s->total)/dnnz;
		if (p > 2.0){
			printf("%2lu: %lu (%.2lf %%) ", (unsigned long)i, s->count, p);
		}
	}
}

static void update_stats(rle_stat_t *stats, unsigned long freq, unsigned long val)
{
	if (freq < FREQ_LIMIT){
		return;
	}

	if (val >= STATS_VAL_NR){
		fprintf(stderr, "darn! I need to put a real hashtable now! :-(\n");
		exit(1);
	}

	stats[val].total += freq;
	stats[val].count ++;
}

static void stats_update_line(rle_stat_t *stats, uint64_t *in, uint64_t insize)
{
	uint64_t rle_freq; // frequency of current value
	uint64_t i, prev, curr;

	if (insize == 0){
		return;
	}

	rle_freq = 1;
	prev = in[0];
	//printf("prev=%lu\n", prev);
	for (i=1; i<insize; i++){
		curr = in[i];
		if (curr == prev){
			rle_freq++;
		} else {
			update_stats(stats, rle_freq, prev);
			rle_freq = 1;
			prev = curr;
		}
	}
	update_stats(stats, rle_freq, prev);
}

static rle_stat_t *rle_stats(point_t *points, uint64_t points_nr)
{
	uint64_t i, prev_y, *xs, *deltas, xs_size, *x;
	rle_stat_t *stats;
	point_t *p;
	dynarray_t *da_xs, *da_deltas;

	da_xs = dynarray_create(sizeof(uint64_t), 512);
	da_deltas = dynarray_create(sizeof(uint64_t), 512);
	stats = calloc(STATS_VAL_NR, sizeof(rle_stat_t));
	if (!stats){
		perror("calloc");
		exit(1);
	}

	void do_update()
	{
		xs_size = dynarray_size(da_xs);
		if (xs_size == 0)
			return;
		xs = dynarray_get(da_xs, 0);
		deltas = dynarray_alloc_nr(da_deltas, xs_size);
		delta_encode(xs, deltas, xs_size);
		stats_update_line(stats, deltas, xs_size);
		dynarray_dealloc_all(da_xs);
		dynarray_dealloc_all(da_deltas);
	}

	prev_y = 1;
	for (i=0; i<points_nr; i++){
		p = points + i;
		// we assume that points are allready sorted
		if (prev_y != p->y){
			do_update();
			prev_y = p->y;
		}
		x = dynarray_alloc(da_xs);
		*x = p->x;
	}
	do_update();

	free(dynarray_destroy(da_xs));
	return stats;
}

static point_t *xform_to_vert(point_t *pts,
                              uint64_t nrows, uint64_t ncols, uint64_t nnz);
{
	point_t *vpts;
	vpts = calloc(nnz, sizeof(point_t));
	if (!vpts){
		perror("calloc");
		exit(1);
	}
	for (i=0; i<nnz; i++){
		vpts[i].x = pts[i].y;
		vpts[i].y = pts[i].x;
	}
	qsort(vpts, nnz, sizeof(point_t), point_cmp_fn);
	return vpts;
}

static point_t *xform_to_diag(point_t *pts,
                              uint64_t nrows, uint64_t ncols, uint64_t nnz);
{
	point_t *dpts;
	dpts = calloc(nnz, sizeof(point_t));
	if (!dpts){
		perror("calloc");
		exit(1);
	}
	for (i=0; i<nnz; i++){
		point_t *p = pts + i;
		uint64_t px = p->x, py=p->y;
		assert(rows + px - py > 0);
		dpts[i].y = nrows + px - py;
		dpts[i].x = (px < py) ? px : py;
	}
	qsort(dpts, nnz, sizeof(point_t), point_cmp_fn);
	return dpts;
}

static point_t *xform_to_rdiag(point_t *pts,
                               uint64_t nrows, uint64_t ncols, uint64_t nnz);
{
	point_t *rdpts;
	rdpts = calloc(nnz, sizeof(point_t));
	if (!rdpts){
		perror("calloc");
		exit(1);
	}
	for (i=0; i<nnz; i++){
		point_t *p = pts + i;
		uint64_t px = p->x, py=p->y;
		uint64_t ny;
		rdpts[i].y = ny = px +py -1;
		rdpts[i].x = ny <= nrows ? px : px +rows -ny;
	}
	qsort(rdpts, nnz, sizeof(point_t), point_cmp_fn);
	return rdpts;
}


spm_t *spm_init_mmf(char *mmf_file)
{
	spm_t *spm;
	spm_elem_t *elem;

	spm = malloc(sizeof(spm_t));
	if (!spm){
		perror("malloc");
		exit(1);
	}

	mmf = mmf_init(mmf_file, &spm->nrows, &spm->ncols, &spm->nnz)
	spm->type = HORIZONTAL;
	spm->elems = malloc(sizeof(spm_elem_t)*nnz);
	if (!spm->elems){
		perror("malloc");
		exit(1);
	}

	elem = spm->elems;
	while (mmf_get_next(mmf, &row, &col, &dummy)){
		elem->type = ELEM_SINGLE;
		elem->point = {.x = row + 1, .y = col + 1};
		elem++;
	}
	spm->elems_nr = elem - spm->elems;
	assert(spm->elems_nr == spm->nnz);
	fclose(mmf);
	return spm;
}

void spm_dealloc(spm_t  *spm)
{
	if (spm->elems)
		free(spm->elems);
	free(spm);
}

spm_t *spm_rle_group(spm_t *spm, long limit)
{
	dynarray_t *elems_da;
	elems_da = dynnarray_create(sizeof(spm_elem_t), 512);
}

int main(int argc, char **argv)
{
}
