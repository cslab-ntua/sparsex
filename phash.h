#ifndef PHASH_H
#define PHASH_H

typedef unsigned long ul_t;
typedef struct {
    ul_t k;
    ul_t v;
} item_t;

struct phash {
    item_t *items;
    ul_t size;
    ul_t used;
    ul_t dummies;
    ul_t minsize;
    ul_t defval;
#ifdef STATS
    ul_t inserts;
    ul_t deletes;
    ul_t lookups;
    ul_t bounces;
#endif
};

#ifdef STATS
#define ZEROSTAT(stat) (stat) = 0
#define INCSTAT(stat) (stat) ++
#define DECSTAT(stat) (stat) --
#else
#define ZEROSTAT(stat)
#define INCSTAT(stat)
#define DECSTAT(stat)
#endif

#define VALID ((~(ul_t)0 >> 1) + 1)
#define DUMMY 1

struct phash *phash_new(ul_t minsize, ul_t defval);
void phash_free(struct phash *phash);
ul_t phash_grow(struct phash *phash);
//ul_t phash_shrink(struct phash *phash);
ul_t phash_insert(struct phash *phash, ul_t key, ul_t val);
//ul_t phash_delete(struct phash *phash, ul_t key);
int phash_lookup(struct phash *phash, ul_t key, ul_t *val);

static inline ul_t phash_elements(struct phash *phash)
{
	return (phash->used - phash->dummies);
}

typedef struct phash phash_t;

#endif
