#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "method.h"

static method_t *methods = NULL ;

method_t *method_get(char *name)
{
	method_t *m;

	for ( m = methods ; m ; m = m->next) {
		if ( strcmp(m->name, name) == 0 )
			break;
	}

	return m;
}

method_t *method_create(char *name, void *fn, void *data)
{
	method_t *ret;
	int name_size = strlen(name);

	ret = malloc( sizeof(*ret) + (name_size+1) );
	if (!ret){
		perror("malloc failed\n");
		return NULL;
	}

	ret->next = NULL;
	ret->fn   = fn;
	ret->data = data;
	ret->name = ((void *)ret + sizeof(*ret));
	strcpy(ret->name, name);

	return ret;
}

void method_add(method_t *method)
{
	method->next = NULL;

	if ( !methods ){
		methods = method;
		return;
	}

	method_t *m;
	for ( m = methods ; m->next ; m = m->next)
		;

	m->next = method;
}

void method_print(char *start, char *sep_start, char *sep_end, char *end)
{
	method_t *m;

	printf("%s", start);
	for ( m = methods ; m ; m = m->next ){
		printf("%s%s%s", sep_start, m->name, sep_end);
	}
	printf("%s", end);
}

void method_fprint(FILE *stream, char *start, char *sep_start, char *sep_end, char *end)
{
	method_t *m;

	fprintf(stream, "%s", start);
	for ( m = methods ; m ; m = m->next ){
		fprintf(stream, "%s%s%s", sep_start, m->name, sep_end);
	}
	fprintf(stream, "%s", end);
}
