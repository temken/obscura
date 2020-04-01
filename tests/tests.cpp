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
#include "DM_Particle.hpp"
#include "DM_Particle_Standard.hpp"
#include "Direct_Detection_Nucleus.hpp"

int main(int argc, char *argv[])
{
	//Starting time
	auto time_start = std::chrono::system_clock::now();
	auto time_start_t = std::chrono::system_clock::to_time_t(time_start);
	std::cout 	<<"Started at " <<std::ctime(&time_start_t)<<std::endl;

	//Import configuration file
	Import_Nuclear_Data();
	Configuration cfg("test.cfg");

	Vector vEarth = Earth_Velocity(0.0);
	std::cout<<"vEarth = "<<In_Units(vEarth,km/sec)<<std::endl;

	//Ending time and computing time
	auto time_end = std::chrono::system_clock::now();
	double durationTotal =1e-6*std::chrono::duration_cast<std::chrono::microseconds>( time_end - time_start ).count();
	std::cout 	<<"\n[Finished in "<< Round(durationTotal,2)<<"s";
	if(durationTotal > 60.0) std::cout <<" ("<<floor(durationTotal/3600.0)<<":"<<floor(fmod(durationTotal/60.0,60.0))<<":"<<floor(fmod(durationTotal,60.0))<<")]."<<std::endl;
	else std::cout <<"]"<<std::endl;
	
	return 0;
}