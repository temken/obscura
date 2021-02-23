#ifndef __DM_Particle_hpp__
#define __DM_Particle_hpp__

#include <random>
#include <string>

#include "Target_Nucleus.hpp"

namespace obscura
{

//1. Base class for a DM particle with virtual functions for the cross sections
class DM_Particle
{
  protected:
	bool low_mass, using_cross_section;

	// Base class implementations
	void Print_Summary_Base(int MPI_rank = 0) const;

	double Sigma_Nucleus_Base(const Isotope& target, double vDM) const;
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

	//Differential cross sections
	virtual double dSigma_dq2_Nucleus(double q, const Isotope& target, double vDM) const { return 0.0; };
	virtual double dSigma_dq2_Electron(double q, double vDM) const { return 0.0; };
	double dSigma_dER_Nucleus(double ER, const Isotope& target, double vDM) const;

	// Reference cross sections
	virtual double Sigma_Proton() const { return 0.0; };
	virtual double Sigma_Neutron() const { return 0.0; };
	virtual double Sigma_Electron() const { return 0.0; };

	virtual double Sigma_Nucleus(const Isotope& target, double vDM) const
	{
		return Sigma_Nucleus_Base(target, vDM);
	};
	virtual double Sigma_Total_Electron(double vDM) const
	{
		return Sigma_Electron_Total_Base(vDM);
	};

	virtual void Print_Summary(int MPI_rank = 0) const
	{
		Print_Summary_Base(MPI_rank);
	};

	// Scattering angle functions
	virtual double PDF_Scattering_Angle_Nucleus(double cos_alpha, const Isotope& target, double vDM)
	{
		return PDF_Scattering_Angle_Nucleus_Base(cos_alpha, target, vDM);
	}

	virtual double PDF_Scattering_Angle_Electron(double cos_alpha, double vDM)
	{
		return PDF_Scattering_Angle_Electron_Base(cos_alpha, vDM);
	}

	virtual double CDF_Scattering_Angle_Nucleus(double cos_alpha, const Isotope& target, double vDM)
	{
		return CDF_Scattering_Angle_Nucleus_Base(cos_alpha, target, vDM);
	}

	virtual double CDF_Scattering_Angle_Electron(double cos_alpha, double vDM)
	{
		return CDF_Scattering_Angle_Electron_Base(cos_alpha, vDM);
	}

	virtual double Sample_Scattering_Angle_Nucleus(const Isotope& target, double vDM, std::mt19937& PRNG)
	{
		return Sample_Scattering_Angle_Nucleus_Base(target, vDM, PRNG);
	}

	virtual double Sample_Scattering_Angle_Electron(double vDM, std::mt19937& PRNG)
	{
		return Sample_Scattering_Angle_Electron_Base(vDM, PRNG);
	}
};

}	// namespace obscura

#endif
