#include <iostream>
#include <chrono>
#include <cmath>

//Headers from libphys library
#include "Natural_Units.hpp"
#include "Linear_Algebra.hpp"

#include "DM_Distribution.hpp"
#include "Astronomy.hpp"
#include "Target_Nucleus.hpp"

int main(int argc, char *argv[])
{
	//Starting time
		std::chrono::high_resolution_clock::time_point time_start = std::chrono::high_resolution_clock::now();

	Standard_Halo_Model SHM;
	SHM.Print_Summary();
	std::cout <<SHM.PDF_Speed(300*km/sec)<<std::endl;
	std::cout <<In_Units(SHM.Average_Speed(),km/sec)<<std::endl;

	Vector vEarth = Earth_Velocity(0.0);
	std::cout<<"vEarth = "<<In_Units(vEarth,km/sec)<<std::endl;

	Import_Nuclear_Data();
	Get_Element(2).Print_Summary();
	
	//Ending time and computing time
	std::chrono::high_resolution_clock::time_point time_end = std::chrono::high_resolution_clock::now();
	double durationTotal =1e-6*std::chrono::duration_cast<std::chrono::microseconds>( time_end - time_start ).count();
	std::cout <<"\nProcessing Time:\t"<< durationTotal<<"s ("<< floor(durationTotal/3600.0)<<":"<<floor(fmod(durationTotal/60.0,60.0))<<":"<<floor(fmod(durationTotal,60.0))<<":"<<floor(fmod(1000*durationTotal,1000.0))<<")."<<std::endl
	<<"##############################"<<std::endl;
	
	return 0;
}