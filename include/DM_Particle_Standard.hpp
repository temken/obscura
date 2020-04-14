#ifndef __DM_Particle_Standard_hpp__
#define __DM_Particle_Standard_hpp__

#include <string>

//Headers from libphys library
#include "Natural_Units.hpp"

#include "DM_Particle.hpp"


//1. Abstract parent class to encompass SI and SD interactions
class DM_Particle_Standard : public DM_Particle
{
	protected:
		//Effective proton and neutron couplings
		double fp, fn;

		//Prefactor of proton/neutron cross section (=1 for SI, and =3 for SD)
		double prefactor;

		//Relation between couplings
		bool fixed_coupling_relation;
		double fp_relative;
		double fn_relative;

		//Reference electron cross section
		double sigma_electron;

		void Print_Summary_Standard(int MPI_rank = 0) const;

	public:
		DM_Particle_Standard();
		DM_Particle_Standard(double mDM,double pre);

		//Primary interaction parameter, in this case the proton or neutron cross section
		virtual double Get_Interaction_Parameter(std::string target) const override;
		virtual void Set_Interaction_Parameter(double par, std::string target) override;

		virtual void Set_Sigma_Proton(double sigma) override;
		virtual void Set_Sigma_Neutron(double sigma) override;
		virtual void Set_Sigma_Electron(double sigma) override;

		void Fix_Coupling_Ratio(double fp_rel, double fn_rel);
		void Fix_fn_over_fp(double ratio);
		void Fix_fp_over_fn(double ratio);
		void Unfix_Coupling_Ratios();
		
		//Differential cross sections with nuclear isotopes, elements, and electrons
		virtual double dSigma_dq2_Nucleus(double q,const Isotope& target,double vDM) const override {return 0;};
		virtual double dSigma_dq2_Electron(double q,double vDM) const override {return 0;};

		//Reference cross sections
		virtual double Sigma_Proton() const override;
		virtual double Sigma_Neutron() const override;
		virtual double Sigma_Electron() const override;

		//Total cross sections with nuclear isotopes, elements, and electrons
		virtual double Sigma_Nucleus(const Isotope& target,double vDM=1e-3) const  override {return 0;};
};

//2. Spin-independent (SI) interactions
class DM_Particle_SI : public DM_Particle_Standard
{
	private:
		double qRef = aEM * mElectron;
		//Dark matter form factor
		std::string FF_DM;
		double mMediator;
		double FormFactor2_DM(double q) const;

	public:
		DM_Particle_SI();
		DM_Particle_SI(double mDM, double sigmaP=1e-40*cm*cm);

		void Set_FormFactor_DM(std::string ff,double mMed=-1.0);
		void Set_Mediator_Mass(double m);
		
		//Differential cross sections with nuclear isotopes, elements, and electrons
		virtual double dSigma_dq2_Nucleus(double q,const Isotope& target,double vDM) const override;
		virtual double dSigma_dq2_Electron(double q,double vDM) const override;

		//Total cross sections with nuclear isotopes, elements, and electrons
		virtual double Sigma_Nucleus(const Isotope& isotope,double vDM=1e-3) const override;

		virtual void Print_Summary(int MPI_rank = 0) const override;
};

//3. Spin-dependent (SD) interactions
class DM_Particle_SD : public DM_Particle_Standard
{
	private:
		//Nuclear form factors missing.

	public:
		DM_Particle_SD();
		DM_Particle_SD(double mDM, double sigmaP=1e-40*cm*cm);
		
		//Differential cross sections with nuclear isotopes, elements, and electrons
		virtual double dSigma_dq2_Nucleus(double q,const Isotope& target,double vDM) const override;
		virtual double dSigma_dq2_Electron(double q,double vDM) const override;

		//Total cross sections with nuclear isotopes, elements, and electrons
		virtual double Sigma_Nucleus(const Isotope& isotope,double vDM=1e-3) const override;

		virtual void Print_Summary(int MPI_rank = 0) const override;
};

#endif