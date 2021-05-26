#ifndef __DM_Particle_hpp__
#define __DM_Particle_hpp__

#include <random>
#include <string>

#include "obscura/Target_Atom.hpp"
#include "obscura/Target_Crystal.hpp"
#include "obscura/Target_Nucleus.hpp"

namespace obscura
{

//1. Base class for a DM particle with virtual functions for the cross sections
class DM_Particle
{
  protected:
	bool low_mass, using_cross_section;

	// Base class implementations
	void Print_Summary_Base(int MPI_rank = 0) const;

	double Sigma_Nucleus_Total_Base(const Isotope& target, double vDM) const;
	double Sigma_Electron_Total_Base(double vDM) const;

	double PDF_Scattering_Angle_Nucleus_Base(double cos_alpha, const Isotope& target, double vDM);
	double PDF_Scattering_Angle_Electron_Base(double cos_alpha, double vDM);
	double CDF_Scattering_Angle_Nucleus_Base(double cos_alpha, const Isotope& target, double vDM);
	double CDF_Scattering_Angle_Electron_Base(double cos_alpha, double vDM);
	double Sample_Scattering_Angle_Nucleus_Base(const Isotope& target, double vDM, std::mt19937& PRNG);
	double Sample_Scattering_Angle_Electron_Base(double vDM, std::mt19937& PRNG);

  public:
	double mass, spin, fractional_density;
	bool DD_use_eta_function;

	//Constructors:
	DM_Particle();
	explicit DM_Particle(double m, double s = 1.0 / 2.0);

	virtual void Set_Mass(double mDM);
	void Set_Spin(double s);
	void Set_Low_Mass_Mode(bool ldm);
	void Set_Fractional_Density(double f);

	//Primary interaction parameter, such as a coupling constant or cross section
	virtual double Get_Interaction_Parameter(std::string target) const
	{
		return 0.0;
	};
	virtual void Set_Interaction_Parameter(double par, std::string target) {};
	bool Interaction_Parameter_Is_Cross_Section() const;

	// Set reference cross sections
	virtual void Set_Sigma_Proton(double sigma) {};
	virtual void Set_Sigma_Neutron(double sigma) {};
	virtual void Set_Sigma_Electron(double sigma) {};

	//Differential cross sections for nuclear targets
	virtual double dSigma_dq2_Nucleus(double q, const Isotope& target, double vDM) const { return 0.0; };
	double dSigma_dER_Nucleus(double ER, const Isotope& target, double vDM) const;
	double d2Sigma_dER_dEe_Migdal(double ER, double Ee, double vDM, const Isotope& isotope, Atomic_Electron& shell) const;

	// Differential cross section for electron targets
	virtual double dSigma_dq2_Electron(double q, double vDM) const { return 0.0; };
	virtual double d2Sigma_dq2_dEe_Ionization(double q, double Ee, double vDM, Atomic_Electron& shell) const { return 0.0; };
	virtual double d2Sigma_dq2_dEe_Crystal(double q, double Ee, double vDM, Crystal& crystal) const { return 0.0; };

	// Reference cross sections
	virtual double Sigma_Proton() const { return 0.0; };
	virtual double Sigma_Neutron() const { return 0.0; };
	virtual double Sigma_Electron() const { return 0.0; };

	virtual double Sigma_Total_Nucleus(const Isotope& target, double vDM) const;
	virtual double Sigma_Total_Electron(double vDM) const;

	virtual void Print_Summary(int MPI_rank = 0) const;

	// Scattering angle functions
	virtual double PDF_Scattering_Angle_Nucleus(double cos_alpha, const Isotope& target, double vDM);
	virtual double PDF_Scattering_Angle_Electron(double cos_alpha, double vDM);
	virtual double CDF_Scattering_Angle_Nucleus(double cos_alpha, const Isotope& target, double vDM);
	virtual double CDF_Scattering_Angle_Electron(double cos_alpha, double vDM);
	virtual double Sample_Scattering_Angle_Nucleus(const Isotope& target, double vDM, std::mt19937& PRNG);
	virtual double Sample_Scattering_Angle_Electron(double vDM, std::mt19937& PRNG);
};

}	// namespace obscura

#endif
