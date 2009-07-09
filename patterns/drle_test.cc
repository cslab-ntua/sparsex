#include <iostream>

#include "spm.h"
#include "drle.h"

using namespace csx;

int main(int argc, char **argv)
{
	SpmIdx *spm;
	DRLE_Manager *drle_mg;
	SpmIterOrder type;

	if (argc < 2){
		std::cerr << "Usage: " << argv[0] << " <mmf_file>" << std::endl;
		return 1;
	}

	spm = loadMMF_mt(argv[1], 1);
	std::cout << "Loaded Matrix " << argv[1] << std::endl;
	std::cout << spm->rows;
	//TestMMF(spm, argv[1]);

	drle_mg = new DRLE_Manager(spm);
	drle_mg->genAllStats();
	drle_mg->outStats();
	std::cout << "Choose Type" << std::endl;
	type = drle_mg->chooseType();
	if (type == NONE)
		return 0;

	std::cout << "Transform to " << SpmTypesNames[type] << std::endl;
	spm->Transform(type);
	std::cout << "Encode" << std::endl;
	drle_mg->Encode();
	std::cout << "Transform to " << SpmTypesNames[HORIZONTAL] << std::endl;
	spm->Transform(HORIZONTAL);
	TestMMF(spm, argv[1]);
	std::cout << spm->rows << std::endl;

	std::cout << "Stats (after encoding)" << std::endl;
	drle_mg = new DRLE_Manager(spm);
	drle_mg->genAllStats();
	drle_mg->outStats();
}
