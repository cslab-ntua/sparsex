#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <unistd.h>

#include "phash.h"

#define PERTURB_SHIFT 5
#define PERTURB_MASK (((ul_t)1 << PERTURB_SHIFT) - 1)

struct phash *phash_new(ul_t minsize, ul_t defval) {
    struct phash *phash;
    phash = malloc(sizeof(struct phash));
    if (!phash) {
        perror("malloc");
        return NULL;
    }
    phash->items = calloc(1, ((ul_t)1 << minsize) * sizeof(item_t));
    if (!phash->items) {
        perror("calloc");
        return NULL;
    }
    phash->used = 0;
    phash->dummies = 0;
    phash->minsize = minsize;
    phash->defval = defval;
    phash->size = minsize;
    ZEROSTAT(phash->inserts);
    ZEROSTAT(phash->deletes);
    ZEROSTAT(phash->lookups);
    ZEROSTAT(phash->bounces);
    return phash;
}

void phash_free(struct phash *phash) {
    free(phash->items);
    free(phash);
}

ul_t phash_grow(struct phash *phash) {
    ul_t old_size = phash->size;
    item_t *new_items, *old_items = phash->items;
    ul_t new_size, u, i;

    u = phash->used;
    if (u/2 + u >= ((ul_t)1 << phash->size)) {
        new_size = old_size + 1;
    } else {
        new_size = old_size;
    }

    #if 0
    if (((ul_t)1 << new_size) > VALID) {
        perror("not enough memory!");
        return old_size;
    }
    #endif

    new_items = calloc(1, ((ul_t)1 << new_size) * sizeof(item_t));
    if (!new_items) {
        perror("calloc");
        return old_size;
    }

    phash->items = new_items;
    phash->used = 0;
    phash->dummies = 0;
    phash->size = new_size;

    //printf("Growing %lu -> %lu\n", old_size, new_size);
    old_size = ((ul_t)1 << old_size);
    for (i = 0; i < old_size; i++) {
        item_t *e;
        e = &old_items[i];
        if (e->v & VALID) {
            phash_insert(phash, e->k, (e->v & ~VALID));
        }
    }
    //printf("Grown %lu -> %lu\n", old_size, (ul_t)1 << new_size);
    free(old_items);

    return new_size;
}

#if 0
ul_t phash_shrink(struct phash *phash) {
    ul_t old_size = phash->size;
    item_t *new_items, *old_items = phash->items;
    ul_t new_size, u, i;

    u = phash->used;
    new_size = old_size - 1;
    if (new_size == 0) {
        new_size = 1;
    }

    if (((ul_t)1 << new_size) > VALID) {
        perror("not enough memory!");
        return old_size;
    }
    new_items = calloc(1, ((ul_t)1 << new_size) * sizeof(item_t));
    if (!new_items) {
        perror("calloc");
        return old_size;
    }

    phash->items = new_items;
    phash->used = 0;
    phash->dummies = 0;
    phash->size = new_size;

    old_size = ((ul_t)1 << old_size);
    for (i = 0; i < old_size; i++) {
        item_t *e;
        e = &old_items[i];
        if (e->k & VALID) {
            phash_insert(phash, e->k, e->v);
        }
    }

    free(old_items);

    return new_size;
}
#endif

ul_t phash_insert(struct phash *phash, ul_t key, ul_t val) {
    register ul_t s = phash->size, u = phash->used + phash->dummies;
    register item_t *items = phash->items;
    register ul_t h, p;

    if (val & VALID){
        fprintf(stderr, "value: %lu has VALID bit set\n", val);
        exit(1);
    }
    register ul_t v = val & ~VALID;

    p = key >> -s;
    s = ((ul_t)1 << s);

    /*
    printf("INSERTING %lx,%lx @ ", key, val);
    ul_t i;
    for (i = 0; i < s; i++) {
        item_t *e;
        e = &phash->items[i];
        printf("%lx:%lx, ", e->k, e->v);
    }
    printf("\n");
    */

    if ((u/2 + u) >= s) {
        s = phash_grow(phash);
        p = key >> -s;
        s = ((ul_t)1 << s);
        items = phash->items;
    }

    s = s - 1;
    h = key & s;

    for (;;) {
        item_t *e;
        ul_t ev;
        e = &items[h];
        ev = e->v;
        if (!(ev & VALID)) {
            e->k = key;
            e->v = v | VALID;
            phash->used ++;
            // XXX phash->dummies -= ek;
            //printf(":insert %lu @ %lu of %lu :: (%lu, %lu)\n", key, h, (ul_t)1 << phash->size, phash->used, phash->dummies);
            INCSTAT(phash->inserts);
            return h;
        }
        //printf("insert (%lx:%lx) bouncing on %lx, (%lx:%lx)\n", key, val, h, e->k, e->v);
        INCSTAT(phash->bounces);
        h = (4*h + h + 1 + p) & s;
        p >>= PERTURB_SHIFT;
    }

    return VALID;
}

#if 0
// removes are not supported
ul_t phash_delete(struct phash *phash, ul_t key) {
    register ul_t s = phash->size, u = phash->used;
    register item_t *items = phash->items;
    register ul_t k = key & ~VALID, h, p;

    p = k >> -s;
    s = ((ul_t)1 << s);

    if ((4*u >= s) && (s > phash->minsize)) {
        s = phash_shrink(phash);
        p = k >> -s;
        s = ((ul_t)1 << s);
        items = phash->items;
    }

    s = s - 1;
    h = k & s;

    for (;;) {
        item_t *e;
        ul_t ek;
        e = &items[h];
        ek = e->k;
        if (ek & VALID) {
            if (ek == (k | VALID)) {
                e->k = DUMMY;
                phash->used --;
                phash->dummies ++;
                INCSTAT(phash->deletes);
                return h;
            }
        } else if (ek != DUMMY) {
            return h;
        }
        INCSTAT(phash->bounces);
        h = (4*h + h + 1 + p) & s;
        p >>= PERTURB_SHIFT;
    }

    return VALID;
}
#endif

int phash_lookup(struct phash *phash, ul_t key, ul_t *val) {

    register ul_t s = phash->size;
    register item_t *items = phash->items;
    register ul_t h, p;

    p = key >> -s;
    s = ((ul_t)1 << s) - 1;
    h = key & s;

    INCSTAT(phash->lookups);
    for (;;) {
        item_t *e;
        ul_t ek, ev;
        e = &items[h];
        ev = e->v;
        ek = e->k;
        //printf("h=%lu k=%lu ek=%lu ev=%lu\n", h, key, ek, ev);
        if (ev & VALID) {
            if (ek == key) {
                *val = (e->v & ~VALID);
                return 1;
            }
        }  else {
            return 0;
        }
        #if 0
         // no removes yet ...
        else if (ek != DUMMY) {
            INCSTAT(phash->lookups);
            return phash->defval;
        }
        #endif
        INCSTAT(phash->bounces);
        h = (4*h + h + 1 + p) & s;
        p >>= PERTURB_SHIFT;
    }
}

#ifdef PHASH_MAIN
#define BUFLEN 1024
int main(int argc, char **argv)
{
    struct phash *ph;
    char *s, buf[BUFLEN], *endptr;
    double dval;
    ul_t val;
    int ret;

    assert(sizeof(double) == sizeof(ul_t));

    ph = phash_new(12, 0);

    for (;;){
        s = fgets(buf, BUFLEN-1, stdin);
        if (s == NULL){
            break;
        }

        dval = strtod(s, &endptr);
        if (*endptr != '\n'){
            fprintf(stderr, "error parsing %s\n", buf);
            exit(1);
        }

        ul_t key = *((ul_t *)&dval);
        //printf("key=0x%lx\n", key);
        ret = phash_lookup(ph, key, &val);
        if (ret){
            //printf("(exists) %lf : %lu\n", dval, val);
            printf("%lu\n", val);
        } else {
            val = phash_elements(ph);
            phash_insert(ph, key, val);
            //printf("(insert) %lf : %lu\n", dval, val);
            printf("%lu\n", val);
        }
    }

    phash_free(ph);
    return 0;
}
#endif

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
