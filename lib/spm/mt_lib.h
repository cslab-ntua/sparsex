#ifndef __MT_LIB_H__
#define __MT_LIB_H__

void setaffinity_oncpu(unsigned int cpu);
void mt_get_options(unsigned int *nr_cpus, unsigned int **cpus);

#endif
