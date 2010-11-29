#include "mmf.h"

#include <fstream>
#include <cstdlib>

using namespace csx;

namespace csx {

void getMmfHeader(const char *mmf_file, uint64_t &nrows, uint64_t &ncols, uint64_t &nnz)
{
	std::ifstream in;

	in.open(mmf_file);
	getMmfHeader(in, nrows, ncols, nnz);
	in.close();
}

void getMmfHeader(std::istream &in, uint64_t &nrows, uint64_t &ncols, uint64_t &nnz)
{
	char buff[512];
	int ret;

	do {
		in.getline(buff, sizeof(buff));
	} while (buff[0] == '#');
	ret = sscanf(buff, "%lu %lu %lu", &nrows, &ncols, &nnz);
	if (ret != 3){
		std::cerr << "mmf header error: sscanf" << std::endl;
		exit(1);
	}
}

} // end csx namespace

MMF::MMF(std::istream &_in) : in(_in)
{

	getMmfHeader(this->in, this->nrows, this->ncols, this->nnz);
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
