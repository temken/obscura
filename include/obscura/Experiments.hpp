#ifndef __Experiments_hpp_
#define __Experiments_hpp_

#include "obscura/Direct_Detection_Crystal.hpp"
#include "obscura/Direct_Detection_ER.hpp"
#include "obscura/Direct_Detection_Migdal.hpp"
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
extern DM_Detector_Ionization_ER XENON10_S2_ER();
extern DM_Detector_Ionization_ER XENON100_S2_ER();
extern DM_Detector_Ionization_ER XENON1T_S2_ER();
extern DM_Detector_Ionization_ER DarkSide50_S2_ER();
//3. Electron recoil experiments - Semiconductor/crystals
extern DM_Detector_Crystal protoSENSEI_at_Surface();
extern DM_Detector_Crystal protoSENSEI_at_MINOS();
extern DM_Detector_Crystal SENSEI_at_MINOS();
extern DM_Detector_Crystal CDMS_HVeV_2018();
extern DM_Detector_Crystal CDMS_HVeV_2020();
//4. Migdal experiments - Ionization
extern DM_Detector_Ionization_Migdal XENON10_S2_Migdal();
extern DM_Detector_Ionization_Migdal XENON100_S2_Migdal();
extern DM_Detector_Ionization_Migdal XENON1T_S2_Migdal();
extern DM_Detector_Ionization_Migdal DarkSide50_S2_Migdal();

}	// namespace obscura

#endif
