#ifndef CSX_CTL_H__
#define CSX_CTL_H__

extern "C" {
	#include "dynarray.h"
}

namespace csx {

#include "ctl_ll.h"

class CtlManager
{
public:
	SpmIdx *spm;

	class PatInfo {
	public:
		uint8_t flag; // flags allocated for this pattern
		uint64_t nr;  // number of non-zero elemenets
		PatInfo(uint8_t flag_, uint64_t nr_): flag(flag_), nr(nr_) {}
		PatInfo(): flag(0), nr(0) {}
	};
	typedef std::map<long,PatInfo> PatMap;

	PatMap patterns;
	uint8_t flag_avail; // current available flag
	bool row_jmps; // does ctl include row_jmps

	// ctl-encoding parse information
	dynarray_t *ctl_da;
	uint64_t last_col;
	bool new_row; // marker of new_row
	uint64_t empty_rows; // number of empty rows since last non-empty row

	CtlManager(SpmIdx *spm_) :
	spm(spm_), flag_avail(0), row_jmps(false), ctl_da(NULL), last_col(0), empty_rows(0) {}

	uint8_t getFlag(long pattern_id, uint64_t nnz);
	uint8_t *mkCtl(uint64_t *ctl_size);
private:
	void doRow(const SpmRowElems &row);
	void updateNewRow(uint8_t *flags);
	void AddXs(std::vector<uint64_t> &xs);
	void AddPattern(const SpmRowElem &elem, uint64_t jmp);
	uint64_t PreparePat(std::vector<uint64_t> &xs, const SpmRowElem &elem);
};

} // end csx namespace

#endif /* CSX_CTL_H__ */
