#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <inttypes.h>
#include <assert.h>

#include "dynarray.h"
#include "mmf.h"

#define STATS_VAL_NR 256
#define FREQ_LIMIT 4

struct point {
	uint64_t x;
	uint64_t y;
};
typedef struct point point_t;

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

static point_t *pts_parse_mmf(FILE *mmf, unsigned long nnz)
{
	point_t *ret;
	unsigned long i=0;
	unsigned long row, col;
	double dummy;

	ret = malloc(sizeof(point_t)*nnz);
	if (!ret){
		perror("malloc");
		exit(1);
	}

	while (mmf_get_next(mmf, &row, &col, &dummy)){
		ret[i].y = row + 1;
		ret[i].x = col + 1;
		i++;
	}

	assert(i == nnz);
	fclose(mmf);
	return ret;
}

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

int main(int argc, char **argv)
{
	unsigned long rows, cols, nnz, i;
	FILE *mmf;
	point_t *pts, *vpts, *dpts, *rdpts;
	rle_stat_t *hstats, *vstats, *dstats, *rdstats;

	if (argc < 2){
		fprintf(stderr, "Usage: %s <mmf_file>\n", argv[0]);
		exit(1);
	}

	/* horizontal */
	mmf = mmf_init(argv[1], &rows, &cols, &nnz);
	pts = pts_parse_mmf(mmf, nnz);
	hstats = rle_stats(pts, nnz);
	/* vertical */
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
	vstats = rle_stats(vpts, nnz);
	free(vpts);
	/* diagonal */
	dpts = calloc(nnz, sizeof(point_t));
	if (!dpts){
		perror("calloc");
		exit(1);
	}
	for (i=0; i<nnz; i++){
		point_t *p = pts + i;
		assert(rows + p->x - p->y > 0);
		dpts[i].y = rows + p->x - p->y;
		dpts[i].x = (p->x < p->y) ? p->x : p->y;
	}
	qsort(dpts, nnz, sizeof(point_t), point_cmp_fn);
	dstats = rle_stats(dpts, nnz);
	free(dpts);
	/* reverse diagonal */
	rdpts = calloc(nnz, sizeof(point_t));
	if (!rdpts){
		perror("calloc");
		exit(1);
	}
	for (i=0; i<nnz; i++){
		point_t *p = pts + i;
		uint64_t ny;
		rdpts[i].y = ny = p->x +p->y -1;
		rdpts[i].x = ny <= rows ? p->x : p->x +rows -ny;
	}
	qsort(rdpts, nnz, sizeof(point_t), point_cmp_fn);
	rdstats = rle_stats(rdpts, nnz);
	free(rdpts);

	printf("%s", argv[1]);
	printf("\nH: ");
	print_stats(hstats, nnz);
	printf("\nV: ");
	print_stats(vstats, nnz);
	printf("\nD: ");
	print_stats(dstats, nnz);
	printf("\nrD: ");
	print_stats(rdstats, nnz);
	printf("\n");

	free(hstats);
	free(vstats);

	return 0;
}
