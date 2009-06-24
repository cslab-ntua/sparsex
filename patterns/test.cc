#include <iostream>
#include <fstream>
#include <cassert>

#include "spm.h"
#include "mmf.h"

using namespace csx;

namespace csx {
void TestMMF(SpmIdx *spm, const char *mmf_file)
{
	SpmIdx::PointIter pi;
	SpmCooElem spm_elem;
	CooElem coo_elem;
	Pattern::Generator *gen;
	std::ifstream in;
	uint64_t y, x;
	double v;

	in.open(mmf_file);
	MMF mmf(in);

	for (pi = spm->points_begin(); pi != spm->points_end(); ++pi){
		spm_elem = *pi;
		coo_elem = static_cast<CooElem>(spm_elem);
		if (spm_elem.pattern == NULL){
			assert(mmf.next(y, x, v));
			assert(coo_elem.y == y);
			assert(coo_elem.x == x);
		} else {
			gen = spm_elem.pattern->generator(coo_elem);
			while ( !gen->isEmpty() ){
				coo_elem = gen->next();
				assert(mmf.next(y, x, v));
				assert(coo_elem.y == y);
				assert(coo_elem.x == x);
			}
		}
	}

	in.close();
}

} // end csx namespace

