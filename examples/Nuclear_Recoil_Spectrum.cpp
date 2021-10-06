#include <iostream>

#include "libphysica/Natural_Units.hpp"

#include "obscura/DM_Halo_Models.hpp"
#include "obscura/DM_Particle_Standard.hpp"
#include "obscura/Direct_Detection_Nucleus.hpp"
#include "obscura/Target_Nucleus.hpp"

using namespace libphysica::natural_units;

int main()
{

	// 1. DM particle (SI and Sd)
	obscura::DM_Particle_SI dm_SI(10.0 * GeV);
	dm_SI.Set_Sigma_Proton(1.0e-40 * cm * cm);
	dm_SI.Print_Summary();

	obscura::DM_Particle_SD dm_SD(10.0 * GeV);
	dm_SD.Set_Sigma_Proton(1.0e-40 * cm * cm);
	dm_SD.Print_Summary();

	// 2. DM distribution
	obscura::Standard_Halo_Model shm;
	shm.Print_Summary();

	// 3. Direct detection targets
	obscura::Nucleus xenon = obscura::Get_Nucleus("Xe");
	xenon.Print_Summary();

	// 4. Evalute the nuclear recoil spectrum
	double E_R		= 1.0 * keV;
	double dRdER_SI = obscura::dRdER_Nucleus(E_R, dm_SI, shm, xenon);
	double dRdER_SD = obscura::dRdER_Nucleus(E_R, dm_SD, shm, xenon);

	std::cout << "SI-interactions: \tdR/dER (1 keV) = " << In_Units(dRdER_SI, 1.0 / kg / year / keV) << " events / kg / year / keV" << std::endl;
	std::cout << "SD-interactions: \tdR/dER (1 keV) = " << In_Units(dRdER_SD, 1.0 / kg / year / keV) << " events / kg / year / keV" << std::endl;

	return 0;
}