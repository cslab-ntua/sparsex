override MAKE_CPPFLAGS += -I./ -I../ -I../../prfcnt
override CPPFLAGS += $(MAKE_CPPFLAGS)
override LDFLAGS += $(MAKE_LDFLAGS)

COMPILE = $(CC) $(CPPFLAGS) $(CFLAGS)
libvbl_deps = spm_vbl.o spm_vbl_mt.o

all: libvbl.o

libvbl.o: $(libvbl_deps)
	$(LD) -i --allow-multiple-definition $(libvbl_deps) -o libvbl.o

spm_vbl.o: spm_vbl.c spm_vbl.h ../spm_crs.h ../macros.h ../spmv_method.h ../vector.h ../mmf.h
	for t in double; do                          \
	  for ci in 32; do                              \
	    $(COMPILE) -DSPM_CRS_BITS=$$ci -DELEM_TYPE=$$t \
	               -o spm_vbl$${ci}_$${t}.o -c $< ;    \
	  done                                             \
	done
	$(LD) -i spm_vbl32_double.o -o spm_vbl.o

spm_vbl_mt.o: spm_vbl_mt.c spm_vbl.h spm_vbl_mt.h ../spm_crs.h ../spm_mt.h ../mt_lib.h ../spm_crs_mt.h ../spmv_method.h
	for t in double; do                             \
	  for ci in 32; do                                 \
	    $(COMPILE) -DSPM_CRS_BITS=$$ci -DELEM_TYPE=$$t    \
	               -o spm_vbl$${ci}_$${t}_mt.o -c $< ;    \
	  done                                                \
	done
	$(LD) -i spm_vbl32_double_mt.o -o spm_vbl_mt.o

clean:
	$(RM) *.o
