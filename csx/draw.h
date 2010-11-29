#ifndef CSX_SPM_DRAW__
#define CSX_SPM_DRAW__

#include "spm.h"

void Draw(csx::SPM &spm,
          const char *filename,
		  int row_start=0, int row_end=0,
		  const int width=600, const int height=600);

#endif /* CSX_SPM_DRAW__ */
