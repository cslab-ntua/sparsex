#ifndef __PHASH_H__
#define __PHASH_H__

typedef unsigned long ul_t;

typedef struct {
    ul_t k;
    ul_t v;
} item_t;

struct phash {
    item_t *items;
    ul_t size_shift;
    ul_t minsize_shift;
    ul_t used;
    ul_t dummies;
    ul_t defval;
#ifdef PHASH_STATS
    ul_t inserts;
    ul_t deletes;
    ul_t lookups;
    ul_t bounces;
#endif
};
typedef struct phash phash_t;

#ifdef PHASH_STATS
#define ZEROSTAT(stat) (stat) = 0
#define INCSTAT(stat) (stat) ++
#define DECSTAT(stat) (stat) --
#define REPSTAT(stat)  fprintf(stderr, "" # stat  " = %lu \n" , stat)
#else
#define ZEROSTAT(stat)
#define INCSTAT(stat)
#define DECSTAT(stat)
#define REPSTAT(stat)  do { } while (0)
#endif

phash_t *phash_new(ul_t minsize_shift);
void phash_free(phash_t *phash);
void phash_insert(phash_t *phash, ul_t key, ul_t val);
int phash_update(phash_t *phash, ul_t key, ul_t val);
void phash_freql_update(phash_t *phash, ul_t key, ul_t val);
int phash_delete(struct phash *phash, ul_t key);
int phash_lookup(phash_t *phash, ul_t key, ul_t *val);
int phash_iterate(struct phash *phash, ul_t *loc, ul_t *key, ul_t *val);
void phash_print(phash_t *phash);

static inline ul_t phash_elements(struct phash *phash)
{
    return phash->used;
}

#endif

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
