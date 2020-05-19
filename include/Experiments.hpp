#ifndef __Experiments_hpp_
#define __Experiments_hpp_

#include "Direct_Detection_Nucleus.hpp"
#include "Direct_Detection_Ionization.hpp"
#include "Direct_Detection_Semiconductor.hpp"

namespace obscura
{
//1. Nuclear recoil experiments
	extern DM_Detector_Nucleus DAMIC_N();
	extern DM_Detector_Nucleus XENON1T_N();
	extern DM_Detector_Nucleus CRESST_II();
	extern DM_Detector_Nucleus CRESST_III();
	extern DM_Detector_Nucleus CRESST_surface();
//2. Electron recoil experiments - Ionization
	extern DM_Detector_Ionization XENON10_e();
	extern DM_Detector_Ionization XENON100_e();
	extern DM_Detector_Ionization XENON1T_e();
	extern DM_Detector_Ionization DarkSide_50_e();
//3. Electron recoil experiments - Semiconductor
	extern DM_Detector_Semiconductor protoSENSEI_at_Surface();
	extern DM_Detector_Semiconductor protoSENSEI_at_MINOS();
	extern DM_Detector_Semiconductor SENSEI_at_MINOS();
	extern DM_Detector_Semiconductor CDMS_HVeV();

}	// namespace obscura

#endif
