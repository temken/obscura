#ifndef __DM_Particle_Standard_hpp__
#define __DM_Particle_Standard_hpp__

#include <string>

#include "libphysica/Natural_Units.hpp"

#include "obscura/DM_Particle.hpp"

namespace obscura
{

// 1. Abstract parent class to encompass SI and SD interactions
class DM_Particle_Standard : public DM_Particle
{
  protected:
	// Effective proton and neutron couplings
	double fp, fn;

	// Prefactor of proton/neutron cross section (=1 for SI, and =3 for SD)
	double prefactor;

	// Relation between couplings
	bool fixed_coupling_relation;
	double fp_relative;
	double fn_relative;

	// Reference electron cross section
	double sigma_electron;

	void Print_Summary_Standard(int MPI_rank = 0) const;

  public:
	DM_Particle_Standard();
	DM_Particle_Standard(double mDM, double pre);

	virtual void Set_Mass(double mDM) override;

	// Primary interaction parameter, in this case the proton, neutron, or electron cross section
	virtual double Get_Interaction_Parameter(std::string target) const override;
	virtual void Set_Interaction_Parameter(double par, std::string target) override;

	virtual void Set_Sigma_Proton(double sigma) override;
	virtual void Set_Sigma_Neutron(double sigma) override;
	virtual void Set_Sigma_Electron(double sigma) override;

	void Fix_Coupling_Ratio(double fp_rel, double fn_rel);
	void Fix_fn_over_fp(double ratio);
	void Fix_fp_over_fn(double ratio);
	void Unfix_Coupling_Ratios();

	// Differential cross sections
	virtual double dSigma_dq2_Nucleus(double q, const Isotope& target, double vDM, double param = -1.0) const override
	{
		return 0;
	};
	virtual double dSigma_dq2_Electron(double q, double vDM, double param = -1.0) const override
	{
		return 0;
	};

	// Reference cross sections
	virtual double Sigma_Proton() const override;
	virtual double Sigma_Neutron() const override;
	virtual double Sigma_Electron() const override;
};

// 2. Spin-independent (SI) interactions
class DM_Particle_SI : public DM_Particle_Standard
{
  private:
	double qRef;
	// Dark matter form factor
	std::string FF_DM;
	double mMediator;
	double FormFactor2_DM(double q) const;

  public:
	DM_Particle_SI();
	explicit DM_Particle_SI(double mDM);
	DM_Particle_SI(double mDM, double sigmaP);

	void Set_FormFactor_DM(std::string ff, double mMed = -1.0);
	void Set_Mediator_Mass(double m);

	//Differential cross sections for nuclear targets
	virtual double dSigma_dq2_Nucleus(double q, const Isotope& target, double vDM, double param = -1.0) const override;

	// Differential cross section for electron targets
	virtual double dSigma_dq2_Electron(double q, double vDM, double param = -1.0) const override;
	virtual double d2Sigma_dq2_dEe_Ionization(double q, double Ee, double vDM, Atomic_Electron& shell) const override;
	virtual double d2Sigma_dq2_dEe_Crystal(double q, double Ee, double vDM, Crystal& crystal) const override;

	// Total cross sections
	virtual double Sigma_Total_Nucleus(const Isotope& isotope, double vDM = 1e-3, double param = -1.0) const override;
	virtual double Sigma_Total_Electron(double vDM, double param = -1.0) const override;

	// Scattering angle functions
	virtual double PDF_Scattering_Angle_Nucleus(double cos_alpha, const Isotope& target, double vDM, double param = -1.0) override;
	virtual double PDF_Scattering_Angle_Electron(double cos_alpha, double vDM, double param = -1.0) override;
	virtual double CDF_Scattering_Angle_Nucleus(double cos_alpha, const Isotope& target, double vDM, double param = -1.0) override;
	virtual double CDF_Scattering_Angle_Electron(double cos_alpha, double vDM, double param = -1.0) override;
	virtual double Sample_Scattering_Angle_Nucleus(std::mt19937& PRNG, const Isotope& target, double vDM, double param = -1.0) override;
	virtual double Sample_Scattering_Angle_Electron(std::mt19937& PRNG, double vDM, double param = -1.0) override;

	virtual void Print_Summary(int MPI_rank = 0) const override;
};

// 3. Spin-dependent (SD) interactions
class DM_Particle_SD : public DM_Particle_Standard
{
  private:
	// Nuclear form factors missing.

  public:
	DM_Particle_SD();
	explicit DM_Particle_SD(double mDM);
	DM_Particle_SD(double mDM, double sigmaP);

	// Differential cross sections with nuclear isotopes, elements, and electrons
	virtual double dSigma_dq2_Nucleus(double q, const Isotope& target, double vDM, double param = -1.0) const override;
	virtual double dSigma_dq2_Electron(double q, double vDM, double param = -1.0) const override;

	// Total cross sections
	virtual double Sigma_Total_Nucleus(const Isotope& isotope, double vDM = 1e-3, double param = -1.0) const override;
	virtual double Sigma_Total_Electron(double vDM, double param = -1.0) const override;

	// Scattering angle functions
	virtual double PDF_Scattering_Angle_Nucleus(double cos_alpha, const Isotope& target, double vDM, double param = -1.0) override;
	virtual double PDF_Scattering_Angle_Electron(double cos_alpha, double vDM, double param = -1.0) override;
	virtual double CDF_Scattering_Angle_Nucleus(double cos_alpha, const Isotope& target, double vDM, double param = -1.0) override;
	virtual double CDF_Scattering_Angle_Electron(double cos_alpha, double vDM, double param = -1.0) override;
	virtual double Sample_Scattering_Angle_Nucleus(std::mt19937& PRNG, const Isotope& target, double vDM, double param = -1.0) override;
	virtual double Sample_Scattering_Angle_Electron(std::mt19937& PRNG, double vDM, double param = -1.0) override;

	virtual void Print_Summary(int MPI_rank = 0) const override;
};

}	// namespace obscura

#endif
