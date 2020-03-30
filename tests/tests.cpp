#include <iostream>
#include <chrono>
#include <cmath>

//Headers from libphys library
#include "Natural_Units.hpp"
#include "Linear_Algebra.hpp"

#include "DM_Distribution.hpp"
#include "Astronomy.hpp"
#include "Target_Nucleus.hpp"
#include "DM_Particle.hpp"
#include "DM_Particle_Standard.hpp"
#include "Direct_Detection_Nucleus.hpp"

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
	Get_Element(54).Print_Summary();

	DM_Particle_SI DM(100.0*GeV, 1.0e-35*cm*cm);
	DM.Fix_fn_over_fp(0.7);
	DM.Print_Summary();
	std::cout<<In_Units(DM.Sigma_Nucleus(Get_Element(54)[5], 300*km/sec),cm*cm)<<std::endl;
	DM.Set_Low_Mass_Mode(true);
	std::cout<<In_Units(DM.Sigma_Nucleus(Get_Element(54)[5], 300*km/sec),cm*cm)<<std::endl;
	
	Detector_Nucleus detector;
	detector.Set_Observed_Signals(100);
	detector.Print_Summary();
	std::cout<<detector.Upper_Bound(DM,SHM)/cm/cm<<std::endl;

	//Ending time and computing time
	std::chrono::high_resolution_clock::time_point time_end = std::chrono::high_resolution_clock::now();
	double durationTotal =1e-6*std::chrono::duration_cast<std::chrono::microseconds>( time_end - time_start ).count();
	std::cout <<"\nProcessing Time:\t"<< durationTotal<<"s ("<< floor(durationTotal/3600.0)<<":"<<floor(fmod(durationTotal/60.0,60.0))<<":"<<floor(fmod(durationTotal,60.0))<<":"<<floor(fmod(1000*durationTotal,1000.0))<<")."<<std::endl
	<<"##############################"<<std::endl;
	
	return 0;
}