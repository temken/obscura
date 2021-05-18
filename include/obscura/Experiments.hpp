#ifndef __Experiments_hpp_
#define __Experiments_hpp_

#include "obscura/Direct_Detection_Ionization.hpp"
#include "obscura/Direct_Detection_Nucleus.hpp"
#include "obscura/Direct_Detection_Semiconductor.hpp"

namespace obscura
{
//1. Nuclear recoil experiments
extern DM_Detector_Nucleus DAMIC_N_2011();
extern DM_Detector_Nucleus XENON1T_N_2017();
extern DM_Detector_Nucleus CRESST_II();
extern DM_Detector_Nucleus CRESST_III();
extern DM_Detector_Nucleus CRESST_surface();
//2. Electron recoil experiments - Ionization
extern DM_Detector_Ionization XENON10_S2();
extern DM_Detector_Ionization XENON100_S2();
extern DM_Detector_Ionization XENON1T_S2();
extern DM_Detector_Ionization DarkSide_50_S2();
//3. Electron recoil experiments - Semiconductor
extern DM_Detector_Semiconductor protoSENSEI_at_Surface();
extern DM_Detector_Semiconductor protoSENSEI_at_MINOS();
extern DM_Detector_Semiconductor SENSEI_at_MINOS();
extern DM_Detector_Semiconductor CDMS_HVeV_2018();
extern DM_Detector_Semiconductor CDMS_HVeV_2020();

}	// namespace obscura

#endif
