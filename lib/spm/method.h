#ifndef __METHOD_H__
#define __METHOD_H__

#include <stdio.h>

struct method_s {
	char            *name;
	void            *fn;
	struct method_s *next;
	void            *data;
};
typedef struct method_s method_t;

method_t *method_get(char *name);
method_t *method_create(char *name, void *fn, void *data);
void method_add(method_t *method);
void method_print(char *start, char *sep_start, char *sep_end, char *end);
void method_fprint(FILE *stream, char *start, char *sep_start, char *sep_end, char *end);

#define METHOD_INIT(func, data) \
void __attribute__((constructor)) func ## _init (void) \
{ \
	method_t *m = method_create( #func, func, data); \
	method_add(m); \
}

#endif /* __METHOD_H__ */
