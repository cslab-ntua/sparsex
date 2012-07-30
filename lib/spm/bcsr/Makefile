override MAKE_CPPFLAGS += -I./ -I../ -I../../prfcnt
override CPPFLAGS += $(MAKE_CPPFLAGS)
override LDFLAGS += $(MAKE_LDFLAGS)

COMPILE = $(CC) $(CPPFLAGS) $(CFLAGS)

mult_deps = mult/generic.o \
	mult/1x2.o \
	mult/1x3.o \
	mult/1x4.o \
	mult/1x5.o \
	mult/1x6.o \
	mult/1x7.o \
	mult/1x8.o \
	mult/2x1.o \
	mult/2x2.o \
	mult/2x3.o \
	mult/2x4.o \
	mult/3x1.o \
	mult/3x2.o \
	mult/3x3.o \
	mult/4x1.o \
	mult/4x2.o \
	mult/5x1.o \
	mult/6x1.o \
	mult/7x1.o \
	mult/8x1.o

libbcsr_deps = spm_bcsr.o spm_bcsr_mt.o blocks.o util.o $(mult_deps)

all: libbcsr.o

libbcsr.o: $(libbcsr_deps)
	$(LD) -i --allow-multiple-definition $(libbcsr_deps) -o libbcsr.o

blocks.o: blocks.c blocks.h spm_bcsr.h ../spm_crs.h util.h
	for t in double; do                          \
	  for ci in 32; do                              \
	    $(COMPILE) -DSPM_CRS_BITS=$$ci -DELEM_TYPE=$$t \
	               -o blocks$${ci}_$${t}.o -c $< ;    \
	  done                                             \
	done
	$(LD) -i blocks32_double.o -o blocks.o

spm_bcsr.o: spm_bcsr.c spm_bcsr.h ../macros.h ../vector.h ../spm_crs.h ../mmf.h ../spmv_method.h spm_bcsr.h blocks.h bitstr.h util.h
	for t in double; do                          \
	  for ci in 32; do                              \
	    $(COMPILE) -DSPM_CRS_BITS=$$ci -DELEM_TYPE=$$t \
	               -o spm_bcsr$${ci}_$${t}.o -c $< ;    \
	  done                                             \
	done
	$(LD) -i spm_bcsr32_double.o -o spm_bcsr.o

spm_bcsr_mt.o: spm_bcsr_mt.c spm_bcsr_mt.h spm_bcsr.h ../macros.h ../vector.h ../spm_crs.h ../mmf.h ../spmv_method.h spm_bcsr.h blocks.h bitstr.h util.h
	for t in double; do                          \
	  for ci in 32; do                              \
	    $(COMPILE) -DSPM_CRS_BITS=$$ci -DELEM_TYPE=$$t \
	               -o spm_bcsr$${ci}_$${t}_mt.o -c $< ;    \
	  done                                             \
	done
	$(LD) -i spm_bcsr32_double_mt.o -o spm_bcsr_mt.o

util.o: util.c util.h
	$(COMPILE) -c $< -o $@

bitstr.h: util.h

mult/%.o: mult/%.c mult/multiply.h
	$(COMPILE) -DSPM_CRS_BITS=32 -DELEM_TYPE=double -o $@ -c $<
	$(LD) -i $@ -o spm_bcsr_multiply.o

clean:
	$(RM) *.o mult/*.o
