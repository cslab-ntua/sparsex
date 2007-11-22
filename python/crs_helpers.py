def crs_size(nr_nz, nr_rows, i_size=4, v_size=8):
	return nr_nz*v_size, (nr_nz+nr_rows)*i_size
