#ifndef __Experiments_hpp_
#define __Experiments_hpp_

#include "obscura/Direct_Detection_Crystal.hpp"
#include "obscura/Direct_Detection_DMe.hpp"
#include "obscura/Direct_Detection_Nucleus.hpp"

namespace obscura
{
//1. Nuclear recoil experiments
extern DM_Detector_Nucleus DAMIC_N_2011();
extern DM_Detector_Nucleus XENON1T_N_2017();
extern DM_Detector_Nucleus CRESST_II();
extern DM_Detector_Nucleus CRESST_III();
extern DM_Detector_Nucleus CRESST_surface();
//2. Electron recoil experiments - Ionization
extern DM_Detector_Ionization_DMe XENON10_S2();
extern DM_Detector_Ionization_DMe XENON100_S2();
extern DM_Detector_Ionization_DMe XENON1T_S2();
extern DM_Detector_Ionization_DMe DarkSide_50_S2();
//3. Electron recoil experiments - Semiconductor/crystals
extern DM_Detector_Crystal protoSENSEI_at_Surface();
extern DM_Detector_Crystal protoSENSEI_at_MINOS();
extern DM_Detector_Crystal SENSEI_at_MINOS();
extern DM_Detector_Crystal CDMS_HVeV_2018();
extern DM_Detector_Crystal CDMS_HVeV_2020();

}	// namespace obscura

#endif
