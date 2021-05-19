#ifndef __Direct_Detection_DMe_hpp_
#define __Direct_Detection_DMe_hpp_

#include "obscura/Direct_Detection_Ionization.hpp"

namespace obscura
{
//1. Event spectra and rates
extern double dRdEe_Ionization_DMe(double Ee, const DM_Particle& DM, DM_Distribution& DM_distr, double m_nucleus, Atomic_Electron& shell);
extern double dRdEe_Ionization_DMe(double Ee, const DM_Particle& DM, DM_Distribution& DM_distr, Atom& atom);

class DM_Detector_Ionization_DMe : public DM_Detector_Ionization
{
  public:
	DM_Detector_Ionization_DMe()
	: DM_Detector_Ionization() {};
	DM_Detector_Ionization_DMe(std::string label, double expo, std::string atom)
	: DM_Detector_Ionization(label, expo, atom) {};
	DM_Detector_Ionization_DMe(std::string label, double expo, std::vector<std::string> atoms, std::vector<double> mass_fractions = {})
	: DM_Detector_Ionization(label, expo, atoms, mass_fractions) {};

	virtual double dRdE_Ionization(double E, const DM_Particle& DM, DM_Distribution& DM_distr, const Nucleus& nucleus, Atomic_Electron& shell) override;
};

};	 // namespace obscura

#endif