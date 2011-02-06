#include <stdio.h>
#include <stdlib.h>

#include "phash.h"

#define LOOPS 1000000

int main(int argc, char **argv) {
    ul_t loops, i, v;
    struct phash *h;

    if (argc > 1) {
        loops = atoi(argv[1]);
    } else {
        loops = LOOPS;
    }

    h = phash_new(2);
    for (i = 10; i < loops; i++) {
	int ret;
        ul_t r;
        //printf("insert(%lx, %lx)\n", i, -i);
        phash_insert(h, i, -i);
        ret = phash_lookup(h, i, &r);
        if (!ret || (r != -i)) {
            printf("insert error in %lx: %lx != %lx\n", i, r, -i);
        }
        //printf(" ->got(%lx, %lx)\n", i, r);
    }
    for (i = 10; i < loops; i++) {
        int ret = phash_lookup(h, i, &v);
        //printf(" ->got(%lu, %lu)\n", i, v);
        if (!ret || (i != -v)) {
            printf("error in %lu: %lu != %lu\n", i, i, -v);
            getchar();
        }
    }

    phash_free(h);

    return 0;
}
