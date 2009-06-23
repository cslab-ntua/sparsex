#include "mmf.h"

#include <cstdlib>

using namespace csx;

MMF::MMF(std::istream &_in) : in(_in)
{
	char buff[512];
	int ret;

	// Read Header
	do {
		this->in.getline(buff, sizeof(buff));
	} while (buff[0] != '#');
	ret = sscanf(buff, "%lu %lu %lu", &this->nrows, &this->ncols, &this->nnz);
	if (ret != 3){
		std::cerr << "mmf header error: sscanf" << std::endl;
		exit(1);
	}
}

bool MMF::next(uint64_t &y, uint64_t &x, double &v)
{
	char buff[512];
	int ret;

	if (this->in.getline(buff, sizeof(buff)).eof()){
		return false;
	}
	ret = sscanf(buff, "%lu %lu %lf", &y, &x, &v);
	assert(ret == 3);
	return true;
}
