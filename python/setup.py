from distutils.core import setup, Extension

import os

dynarray_dir = os.popen('rsrc resource dynarray').read()[:-1]
dynarray_dep = dynarray_dir + '/dynarray.o'

setup(
	name="pyspm", 
	version="0.001", 
	ext_modules=[ Extension(
		"pyspm", 
		["pyspm.c"],
		include_dirs=["../", dynarray_dir],
		extra_objects = ["../libspmv.o", dynarray_dep],
		extra_link_args = ['-Xlinker', '--allow-multiple-definition'],

	)]
)
