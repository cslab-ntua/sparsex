#include <iostream>
#include <string>

#include <cstdlib>

#include <boost/algorithm/string.hpp>

#include "spm.h"
#include "drle.h"
#include "draw.h"

using namespace csx;

int main(int argc, char **argv)
{

	SpmIdx *spm;
	DRLE_Manager *drle_mg;
	DeltaRLE::Stats drle_stats;

	if (argc < 2){
		std::cerr << "Usage: " << argv[0] << " <mmf_file>\n";
		return 1;
	}

	spm = loadMMF_mt(argv[1], 1);
	std::cout << "==> Loaded Matrix" << argv[1] << std::endl;
	drle_mg = new DRLE_Manager(spm);
	drle_mg->genAllStats();
	drle_mg->outStats();

	#if 0
	std::string input;
	std::vector <std::string> tokens;

	char buff[1024];
	for (;;){
		std::cout << argv[1] << "> ";
		std::cin.getline(buff, sizeof(buff));
		if (std::cin.eof())
			break;

		boost::split(tokens, buff, boost::is_any_of(" "));
		if (tokens[0] == "draw"){
			long row_s, row_e;
			row_s = (tokens.size() > 1) ? atol(tokens[1].c_str()) : 0;
			row_e = (tokens.size() > 2) ? atol(tokens[2].c_str()) : 0;
			std::cout << "Drawing (" << row_s << "->" << row_e << ")\n";
			Draw(drle_mg.spm, "test.png", row_s, row_e);
			system("sh -c 'eog test.png &'");
		} else if (tokens[0] == "xform"){
		} else if (tokens[0] == "stats"){
			;
		} else {
			std::cout << "Unknown command "  << tokens[0] << "\n";
		}
	}
	//drle_stats = drle_mg.generateStats(spm);
	//DRLE_OutStats(drle_stats, spm, std::cout);
	//std::cout << std::endl;
	#endif
	delete[] spm;
	return 0;
}
