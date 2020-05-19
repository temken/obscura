#include "DM_Particle.hpp"

#include <iostream>
#include <cmath>
#include <functional>

//Headers from libphysica library
#include "Natural_Units.hpp"
#include "Numerics.hpp"

namespace obscura
{
	using namespace libphysica::natural_units;

//1. Base class for a DM particle with virtual functions for the cross sections
	DM_Particle::DM_Particle()
	: low_mass(false), mass(10.0*GeV), spin(1.0/2.0), fractional_density(1.0)
	{
	}
	
	DM_Particle::DM_Particle(double m, double s)
	: low_mass(false), mass(m), spin(s), fractional_density(1.0)
	{
	}
	
	void DM_Particle::Set_Mass(double mDM)
	{
		mass = mDM;
	}

	void DM_Particle::Set_Spin(double s)
	{
		spin = s;
	}

	void DM_Particle::Set_Low_Mass_Mode(bool ldm)
	{
		low_mass = ldm;
	}

	void DM_Particle::Set_Fractional_Density(double f)
	{
		fractional_density = f;
	}


	double DM_Particle::Sigma_Nucleus_Base(const Isotope& target, double vDM) const
	{
		//Numerically integrate the differential cross section
		double q2min = 0;
		double q2max = 4.0 * pow(libphysica::Reduced_Mass(mass,target.mass)*vDM,2.0);
		std::function<double(double)> dodq2 = [this,&target,vDM] (double q2)
		{
			return dSigma_dq2_Nucleus(sqrt(q2),target,vDM);
		};
	  	double eps = libphysica::Find_Epsilon(dodq2,q2min,q2max,1.0e-6);
	  	double sigmatot = libphysica::Integrate(dodq2,q2min,q2max,eps);
		return sigmatot;
	}

	void DM_Particle::Print_Summary_Base(int MPI_rank) const
	{
		if(MPI_rank == 0)
		{
			std::cout 	<<std::endl
			<<"----------------------------------------"<<std::endl
			<<"DM particle summary:"<<std::endl;

			double massunit = (mass<keV)? eV: ( (mass<MeV)? keV : ((mass<GeV)? MeV : GeV) );
			std::string massunitstr = (mass<keV)? "eV": ( (mass<MeV)? "keV" : ((mass<GeV)? "MeV" : "GeV") );
			std::cout 	<<"\tMass:\t\t\t" <<In_Units(mass,massunit)<<" "<<massunitstr<<std::endl
						<<"\tSpin:\t\t\t" <<spin<<std::endl
						<<"\tLow mass:\t\t" <<((low_mass)? "[x]" : "[ ]")	<<std::endl;
		}
	}

	double DM_Particle::dSigma_dER_Nucleus(double ER,const Isotope& target,double vDM) const
	{
		double q = sqrt(2.0*target.mass*ER);
		return 2.0 * target.mass * dSigma_dq2_Nucleus(q,target,vDM);
	}

}	// namespace obscura