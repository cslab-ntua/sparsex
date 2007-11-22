from distutils.core import setup, Extension

setup(
	name="pyspm", 
	version="0.001", 
	ext_modules=[ Extension(
		"pyspm", 
		["pyspm.c"],
		include_dirs=["../"],
		extra_objects = ["../libspmv.o", "../dynarray.o"]
	)]
)
