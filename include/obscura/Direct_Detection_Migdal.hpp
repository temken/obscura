#ifndef __Direct_Detection_Migdal_hpp_
#define __Direct_Detection_Migdal_hpp_

#include "obscura/Direct_Detection_Ionization.hpp"

namespace obscura
{
//1. Event spectra and rates
extern double dRdEe_Ionization_Migdal(double Ee, const DM_Particle& DM, DM_Distribution& DM_distr, const Isotope& isotope, Atomic_Electron& shell);
extern double dRdEe_Ionization_Migdal(double Ee, const DM_Particle& DM, DM_Distribution& DM_distr, const Nucleus& nucleus, Atomic_Electron& shell);
extern double dRdEe_Ionization_Migdal(double Ee, const DM_Particle& DM, DM_Distribution& DM_distr, Atom& atom);

//2. Detector class for ionization experiments from DM-electron scatterings.
class DM_Detector_Ionization_Migdal : public DM_Detector_Ionization
{
  public:
	DM_Detector_Ionization_Migdal();
	DM_Detector_Ionization_Migdal(std::string label, double expo, std::string atom);
	DM_Detector_Ionization_Migdal(std::string label, double expo, std::vector<std::string> atoms, std::vector<double> mass_fractions = {});

	virtual double dRdE_Ionization(double E, const DM_Particle& DM, DM_Distribution& DM_distr, const Nucleus& nucleus, Atomic_Electron& shell) override;
};

}	// namespace obscura

#endif