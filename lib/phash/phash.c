#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <unistd.h>

#include "phash.h"

#define UNUSED (~(ul_t)0)      /* this entry was never used */
#define DUMMY  ((~(ul_t)0)-1)  /* this entry was used, but now its empty */

//#define VAL_OVERLOAD
//#define KEY_OVERLOAD
//#define NO_OVERLOAD /* use separate bitarray -- not implemented */

#define PERTURB_SHIFT 5

item_t *items_new(ul_t nr_items)
{
    ul_t i;
    item_t *items = malloc(nr_items*sizeof(item_t));

    if (!items) {
        perror("malloc");
        exit(1);
    }

    for (i=0; i < nr_items; i++){
        #ifdef VAL_OVERLOAD
        items[i].v = UNUSED;
        #elif  KEY_OVERLOAD
        items[i].k = UNUSED;
        #endif
    }

    return items;
}

static inline unsigned item_dummy(item_t *item)
{
    #ifdef VAL_OVERLOAD
    return (item->v == DUMMY);
    #elif KEY_OVERLOAD
    return (item->k == DUMMY);
    #endif
}

static void make_item_dummy(item_t *item)
{
    #ifdef VAL_OVERLOAD
    item->v = DUMMY;
    #elif KEY_OVERLOAD
    item->k = DUMMY;
    #endif
}

static inline unsigned item_unused(item_t *item)
{
    #ifdef VAL_OVERLOAD
    return (item->v == UNUSED);
    #elif KEY_OVERLOAD
    return (item->k == UNUSED);
    #endif
}

static inline unsigned item_valid(item_t *item)
{
    return !(item_dummy(item) || item_unused(item));
}

static inline void assert_key(ul_t key)
{
    #if KEY_OVERLOAD
    assert((key != UNUSED) && (key != DUMMY));
    #endif
}

static inline void assert_val(ul_t val)
{
    #ifdef VAL_OVERLOAD
    assert((val != UNUSED) && (val != DUMMY));
    #endif
}

struct phash *phash_new(ul_t minsize_shift) {
    struct phash *phash;

    phash = malloc(sizeof(struct phash));
    if (!phash) {
        perror("malloc");
        exit(1);
    }

    phash->items = items_new(1UL<<minsize_shift);
    phash->dummies = phash->used = 0;
    phash->size_shift = phash->minsize_shift = minsize_shift;

    ZEROSTAT(phash->inserts);
    ZEROSTAT(phash->deletes);
    ZEROSTAT(phash->lookups);
    ZEROSTAT(phash->bounces);

    return phash;
}

void phash_free(struct phash *phash)
{

    REPSTAT(phash->inserts);
    REPSTAT(phash->deletes);
    REPSTAT(phash->lookups);
    REPSTAT(phash->bounces);

    free(phash->items);
    free(phash);
}


ul_t phash_resize(struct phash *phash, ul_t new_size_shift)
{
    ul_t old_size_shift = phash->size_shift;
    ul_t new_size = (ul_t)1<<new_size_shift;
    item_t *new_items, *old_items;
    //fprintf(stderr, "resize: %lu -> %lu\n", old_size_shift, new_size_shift);

    old_items = phash->items;
    new_items = phash->items = items_new(new_size);
    phash->dummies = phash->used = 0;
    phash->size_shift = new_size_shift;

    ul_t old_size = ((ul_t)1 << old_size_shift);
    ul_t i;
    for (i = 0; i < old_size; i++) {
        item_t *item = &old_items[i];
        if (item_valid(item)){
            //fprintf(stderr, "rs: inserting (%lu,%lu)\n", item->k, item->v);
            phash_insert(phash, item->k, item->v);
        }
    }

    free(old_items);
    return new_size_shift;
}

ul_t phash_grow(struct phash *phash)
{
    ul_t old_size_shift = phash->size_shift;
    ul_t new_size_shift;
    ul_t u;

    u = phash->used;
    if (u/2 + u >= ((ul_t)1 << old_size_shift)) {
        new_size_shift = old_size_shift + 1;
    } else {
        new_size_shift = old_size_shift;
    }

    return phash_resize(phash, new_size_shift);
}

ul_t phash_shrink(struct phash *phash)
{
    ul_t old_size_shift = phash->size_shift;
    ul_t new_size_shift;

    new_size_shift = old_size_shift - 1;
    if (new_size_shift < phash->minsize_shift) {
        new_size_shift = phash->minsize_shift;
    }

    return phash_resize(phash, new_size_shift);
}

static inline ul_t phash_grow_check(phash_t *phash)
{
    ul_t size_shift = phash->size_shift;
    ul_t u = phash->used + phash->dummies;
    ul_t size = (ul_t)1<<size_shift;
    if ((u/2 + u) >= size) {
        size_shift = phash_grow(phash);
        size = (ul_t)1<<size_shift;
    }

    return size;
}

#define PHASH_UPDATE(phash, key, val)               \
{                                                   \
    ul_t size = 1<<(phash->size_shift);             \
    item_t *items = phash->items;                   \
    ul_t perturb = key;                             \
    ul_t mask = size-1;                             \
    ul_t i = key & mask;                            \
                                                    \
    INCSTAT(phash->inserts);                        \
    for (;;) {                                      \
        item_t *item = &items[i];                   \
        if ( !item_valid(item) ){                   \
             PHUPD_SET__(phash, item, key, val);    \
             break;                                 \
        }                                           \
        if (item->k == key){                        \
            PHUPD_UPDATE__(phash, item, key, val);  \
            break;                                  \
        }                                           \
                                                    \
        INCSTAT(phash->bounces);                    \
        i = ((i<<2) + i + 1 + perturb) & mask;      \
        perturb >>= PERTURB_SHIFT;                  \
    }                                               \
}

void static inline phash_upd_set(phash_t *p, item_t *i, ul_t key, ul_t val)
{
    if (item_dummy(i))
        p->dummies--;
    p->used++;
    i->v = val;
    i->k = key;
}

void phash_insert(struct phash *phash, ul_t key, ul_t val)
{

    //fprintf(stderr, "insert: (%lu,%lu)\n", key, val);
    assert_key(key);
    assert_val(val);
    phash_grow_check(phash);
    #define PHUPD_UPDATE__(_p, _i, _k, _v) _i->v = _v
    #define PHUPD_SET__(_p, _i, _k, _v) phash_upd_set(_p, _i, _k, _v)
    PHASH_UPDATE(phash, key, val)
    #undef PHUPD_UPDATE__
    #undef PHUPD_SET__
}


void phash_freql_update(struct phash *phash, ul_t key, ul_t val)
{
    assert_key(key);
    assert_val(val);
    phash_grow_check(phash);
    #define PHUPD_UPDATE__(_p, _i, _k, _v) _i->v += _v
    #define PHUPD_SET__(_p, _i, _k, _v) phash_upd_set(_p, _i, _k, _v)
    PHASH_UPDATE(phash, key, val)
    #undef PHUPD_UPDATE__
    #undef PHUPD_SET__
}

/*
 * note that his function does not modify the internal structure of the hash
 * and thus its safe to use it for updating values during a phash_iterate()
 */
int phash_update(struct phash *phash, ul_t key, ul_t val) {

    //fprintf(stderr, "update: (%lu,%lu)\n", key, val);
    assert_key(key);
    assert_val(val);

    #define PHUPD_UPDATE__(_p, _i, _k, _v) _i->v = _v
    #define PHUPD_SET__(_p, _i, _k, _v) return 0;
    PHASH_UPDATE(phash, key, val)
    #undef PHUPD_UPDATE__
    #undef PHUPD_SET__

    return 1;
}

int phash_delete(struct phash *phash, ul_t key)
{
    assert_key(key);

    ul_t size_shift = phash->size_shift;
    ul_t size = (ul_t)1<<size_shift;

    ul_t u = phash->used;
    if (4*u < size){
        size_shift = phash_shrink(phash);
        size = (ul_t)1<<size_shift;
    }

    item_t *items = phash->items;
    ul_t perturb = key;
    ul_t mask = size-1;
    ul_t i = key & mask;

    for (;;) {
        item_t *item = &items[i];

        if ( item_unused(item) ){
            assert(0);
            return 0;
        }

        if ( !item_dummy(item) && item->k == key){
            INCSTAT(phash->deletes);
            make_item_dummy(item);
            phash->dummies++;
            //fprintf(stderr, "rm: used: %lu\n", phash->used);
            phash->used--;
            return 1;
        }

        INCSTAT(phash->bounces);
        i = ((i<<2) + i + 1 + perturb) & mask;
        perturb >>= PERTURB_SHIFT;
    }
}

int phash_lookup(struct phash *phash, ul_t key, ul_t *val)
{
    assert_key(key);

    ul_t size_shift = phash->size_shift;
    ul_t size = (ul_t)1<<size_shift;
    item_t *items = phash->items;
    ul_t perturb = key;
    ul_t mask = size-1;
    ul_t i = key & mask;

    INCSTAT(phash->lookups);
    for (;;) {
        item_t *item = &items[i];

        if ( item_unused(item) )
            return 0;

        if ( !item_dummy(item) && item->k == key){
            *val = item->v;
            return 1;
        }

        INCSTAT(phash->bounces);
        i = ((i<<2) + i + 1 + perturb) & mask;
        perturb >>= PERTURB_SHIFT;
    }
}

int phash_iterate(struct phash *phash, ul_t *loc, ul_t *key, ul_t *val)
{
    ul_t i = *loc;
    ul_t size = (ul_t)1<<phash->size_shift;
    item_t *items = phash->items;
    INCSTAT(phash->lookups);
    for (;;){
        if (i >= size)
            return 0;
        item_t *item = &items[i++];
        if (item_valid(item)){
            *loc = i;
            *key = item->k;
            *val = item->v;
            return 1;
        }
    }
}

void phash_print(phash_t *phash)
{
    ul_t i, key, val;
    int ret;

    printf("PHASH(%p):\n", phash);
    for (i=0 ;;){
        ret = phash_iterate(phash, &i, &key, &val);
        if (!ret){
            break;
        }
        printf(" 0x%017lx : 0x%017lx\n", key, val);
    }
    printf("\n");
}

#ifdef PHASH_MAIN
#define BUFLEN 1024
void help()
{
    printf("Help:\n"
           "  insert : I <key> <val> \n"
           "  update : U <key> <val> (->v += val if exists) \n"
           "  get    : G <key>       \n"
           "  delete : D <key>       \n"
           "  size   : S             \n"
           "  print  : P             \n");
}

int main(int argc, char **argv)
{
    struct phash *ph;
    char *s, buf[BUFLEN];
    ul_t key, val;
    int ret;

    ph = phash_new(2);

    for (;;){
        s = fgets(buf, BUFLEN-1, stdin);
        if (s == NULL){
            break;
        }

        switch (*s) {
            case 'I':
            ret = sscanf(s+1, "%lu %lu", &key, &val);
            if (ret == 2){
                phash_insert(ph, key, val);
            }
            break;

            case 'U':
            ret = sscanf(s+1, "%lu %lu", &key, &val);
            if (ret == 2){
                phash_freql_update(ph, key, val);
            }
            break;

            case 'G':
            ret = sscanf(s+1, "%lu", &key);
            if (ret == 1){
                ret = phash_lookup(ph, key, &val);
                if (ret){
                    printf("%lu\n", val);
                } else {
                    printf("<None>\n");
                }
            }
            break;

            case 'D':
            ret = sscanf(s+1, "%lu", &key);
            if (ret == 1){
                phash_delete(ph, key);
            }
            break;

            case 'S':
            printf("%lu\n", phash_elements(ph));
            break;

            case 'P':
            phash_print(ph);
            break;

            case '#':
            break;

            default:
            help();
            break;

        }
        fflush(stdout);
    }

    phash_free(ph);
    return 0;
}
#endif
// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
