#ifndef __BEANBAG_H__
#define __BEANBAG_H__

#warning "beanbags are deprecated, use dynarrays"

struct beanbag;
typedef struct beanbag beanbag_t;

void *beanbag_create(unsigned long bean_size, unsigned long beans_grain);
void *beanbag_alloc(struct beanbag *bb);
void *beanbag_alloc_nr(struct beanbag *bb, unsigned long nr);
void *beanbag_get(struct beanbag *bb, unsigned long idx);
void *beanbag_get_last(struct beanbag *bb);
unsigned long beanbag_size(struct beanbag *bb);
void *beanbag_destroy(struct beanbag *bb);

#endif
