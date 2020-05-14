#include <iostream>
#include <chrono>
#include <cmath>

//Headers from libphysica library
#include "Natural_Units.hpp"
#include "Numerics.hpp"
#include "Utilities.hpp"

#include "Configuration.hpp"

int main(int argc, char *argv[])
{
	//Starting time
	auto time_start = std::chrono::system_clock::now();
	auto time_start_t = std::chrono::system_clock::to_time_t(time_start);
	std::cout 	<<"Started at " <<std::ctime(&time_start_t)<<std::endl;

	//Import nuclear data and configuration file
	Import_Nuclear_Data();
	Configuration cfg(argv[1]);
	cfg.Print_Summary();

	std::vector<double> DM_masses = Log_Space(cfg.constraints_mass_min, cfg.constraints_mass_max, cfg.constraints_masses);
	std::vector<std::vector<double>> exclusion_limits = cfg.DM_detector->Upper_Limit_Curve(*(cfg.DM), *(cfg.DM_distr), DM_masses, cfg.constraints_certainty);
	Export_Table(TOP_LEVEL_DIR "results/" + cfg.ID + "/constraints.txt", exclusion_limits,{GeV,cm*cm});
	
	//Ending time and computing time
	auto time_end = std::chrono::system_clock::now();
	double durationTotal =1e-6*std::chrono::duration_cast<std::chrono::microseconds>( time_end - time_start ).count();
	std::cout 	<<"\n[Finished in "<< Round(durationTotal,2)<<"s";
	if(durationTotal > 60.0) std::cout <<" ("<<floor(durationTotal/3600.0)<<":"<<floor(fmod(durationTotal/60.0,60.0))<<":"<<floor(fmod(durationTotal,60.0))<<")]."<<std::endl;
	else std::cout <<"]"<<std::endl;
	
	return 0;
}