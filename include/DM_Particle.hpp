#ifndef __DM_Particle_hpp__
#define __DM_Particle_hpp__

#include "Target_Nucleus.hpp"

//1. Base class for a DM particle with virtual functions for the cross sections
class DM_Particle
{
	protected:
		bool low_mass;

		void Print_Summary_Base(int MPI_rank = 0) const;
		double Sigma_Nucleus_Base(const Isotope& target, double vDM) const;

	public:
		double mass, spin, fractional_density;

		//Constructors:
		DM_Particle();
		DM_Particle(double m, double s = 1.0/2.0);

		void Set_Mass(double mDM);
		void Set_Spin(double s);
		void Set_Low_Mass_Mode(bool ldm);
		void Set_Fractional_Density(double f);

		//Primary interaction parameter, such as a coupling constant or cross section
		virtual double Get_Interaction_Parameter() const {return 0.0;};
		virtual void Set_Interaction_Parameter(double par) {};
		
		virtual void Set_Sigma_Proton(double sigma) {};
		virtual void Set_Sigma_Neutron(double sigma) {};
		virtual void Set_Sigma_Electron(double sigma) {};

		//Differential cross sections
		virtual double dSigma_dq2_Nucleus(double q, const Isotope& target, double vDM) const {return 0.0;};
		virtual double dSigma_dq2_Electron(double q, double vDM) const {return 0.0;};
		double dSigma_dER_Nucleus(double ER,const Isotope& target,double vDM) const;
		
		//Cross sections
		virtual double Sigma_Proton() const {return 0.0;};
		virtual double Sigma_Neutron() const {return 0.0;};
		virtual double Sigma_Electron() const {return 0.0;};
		
		virtual double Sigma_Nucleus(const Isotope& target, double vDM) const {return Sigma_Nucleus_Base(target,vDM);};

		virtual void Print_Summary(int MPI_rank = 0) const {Print_Summary_Base(MPI_rank);};
};

#endif