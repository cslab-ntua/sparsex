.PHONY: all

all: pyspm.so

pyspm.so: pyspm.c ../libspmv.o
	make clean
	python2.5 setup.py install_lib -d .

clean:
	rm -rf pyspm.so build *.pyc
