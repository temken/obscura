#include <iostream>
#include <chrono>
#include <cmath>

//Headers from libphys library
#include "Natural_Units.hpp"
#include "Numerics.hpp"
#include "Linear_Algebra.hpp"

#include "DM_Distribution.hpp"
#include "Astronomy.hpp"
#include "Configuration.hpp"
#include "Target_Nucleus.hpp"
#include "Target_Electron.hpp"
#include "DM_Particle.hpp"
#include "DM_Particle_Standard.hpp"
#include "Direct_Detection_Nucleus.hpp"
#include "Direct_Detection_Semiconductor.hpp"

int main(int argc, char *argv[])
{
	//Starting time
	auto time_start = std::chrono::system_clock::now();
	auto time_start_t = std::chrono::system_clock::to_time_t(time_start);
	std::cout 	<<"Started at " <<std::ctime(&time_start_t)<<std::endl;

	//Import nuclear data and configuration file
	Import_Nuclear_Data();
	Configuration cfg("test.cfg");
	cfg.Print_Summary();
	
	Vector vEarth = Earth_Velocity(0.0);
	std::cout<<"vEarth = "<<In_Units(vEarth,km/sec)<<std::endl<<std::endl;

	std::vector<std::vector<double>> exclusion_limits = cfg.DM_detector->Upper_Limit_Curve(*(cfg.DM), *(cfg.DM_distr), cfg.constraints_mass_min, cfg.constraints_mass_max, cfg.constraints_masses, cfg.constraints_certainty);
	for(int i = 0; i < exclusion_limits.size(); i++)
	{
		std::cout <<i+1 <<")\t" <<Round(exclusion_limits[i][0]) <<" GeV\t" <<Round(In_Units(exclusion_limits[i][1],cm*cm)) <<" cm^2" <<std::endl;
		cfg.DM->Set_Mass(exclusion_limits[i][0]);
		cfg.DM->Set_Interaction_Parameter(exclusion_limits[i][1],"Nuclei");
		std::cout<<"N = "<<cfg.DM_detector->DM_Signals_Total(*(cfg.DM),*(cfg.DM_distr))<<std::endl;
	}
	

	//Ending time and computing time
	auto time_end = std::chrono::system_clock::now();
	double durationTotal =1e-6*std::chrono::duration_cast<std::chrono::microseconds>( time_end - time_start ).count();
	std::cout 	<<"\n[Finished in "<< Round(durationTotal,2)<<"s";
	if(durationTotal > 60.0) std::cout <<" ("<<floor(durationTotal/3600.0)<<":"<<floor(fmod(durationTotal/60.0,60.0))<<":"<<floor(fmod(durationTotal,60.0))<<")]."<<std::endl;
	else std::cout <<"]"<<std::endl;
	
	return 0;
}