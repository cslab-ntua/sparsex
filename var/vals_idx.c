#include <Python.h>

#include <stdio.h>
#include <mcheck.h>

#define BUFLEN 128

int main(int argc, char **argv)
{
	long idx;
	char buf[BUFLEN], *s;

	mtrace();

	PyObject *dict;
	if ( !(dict = PyDict_New()) ){
		fprintf(stderr, "cant create python dict\n");
		exit(1);
	}

	for (;;){
		PyObject *pystr, *pyidx;

		s = fgets(buf, BUFLEN-1, stdin);
		if (s == NULL){
			break;
		}

		s[strlen(s) - 1] = '\0'; // remove newline
		pystr = PyString_FromString(s);

		pyidx = PyDict_GetItem(dict, pystr);
		if ( pyidx == NULL){
			idx = PyDict_Size(dict);
			pyidx = PyLong_FromLong(idx);
			PyDict_SetItem(dict, pystr, pyidx);
			Py_DECREF(pystr);
			Py_DECREF(pyidx);
		} else {
			idx = PyLong_AsLong(pyidx);
			Py_DECREF(pystr);
		}

		printf("%ld %s\n", idx, s);
	}

	printf("destroying dict ...\n");
	PyDict_Clear(dict);
	Py_DECREF(dict);

	sleep(1);


	return 0;
}
