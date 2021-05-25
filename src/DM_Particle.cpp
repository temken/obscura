#include "obscura/DM_Particle.hpp"

#include <cmath>
#include <functional>
#include <iostream>

#include "libphysica/Natural_Units.hpp"
#include "libphysica/Numerics.hpp"
#include "libphysica/Statistics.hpp"

namespace obscura
{
using namespace libphysica::natural_units;

//1. Base class for a DM particle with virtual functions for the cross sections
DM_Particle::DM_Particle()
: low_mass(false), using_cross_section(false), mass(10.0 * GeV), spin(1.0 / 2.0), fractional_density(1.0), DD_use_eta_function(false)
{
}

DM_Particle::DM_Particle(double m, double s)
: low_mass(false), using_cross_section(false), mass(m), spin(s), fractional_density(1.0), DD_use_eta_function(false)
{
}

void DM_Particle::Set_Mass(double mDM)
{
	// The cross sections might depend on the DM mass, and need to be re-computed to yield the same cross sections.
	double sigma_p = Sigma_Proton();
	double sigma_n = Sigma_Neutron();
	double sigma_e = Sigma_Electron();

	mass = mDM;

	Set_Sigma_Proton(sigma_p);
	Set_Sigma_Neutron(sigma_n);
	Set_Sigma_Electron(sigma_e);
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

bool DM_Particle::Interaction_Parameter_Is_Cross_Section() const
{
	return using_cross_section;
}

double DM_Particle::Sigma_Nucleus_Total_Base(const Isotope& target, double vDM) const
{
	//Numerically integrate the differential cross section
	double q2min						= 0;
	double q2max						= 4.0 * pow(libphysica::Reduced_Mass(mass, target.mass) * vDM, 2.0);
	std::function<double(double)> dodq2 = [this, &target, vDM](double q2) {
		return dSigma_dq2_Nucleus(sqrt(q2), target, vDM);
	};
	double eps		= libphysica::Find_Epsilon(dodq2, q2min, q2max, 1.0e-6);
	double sigmatot = libphysica::Integrate(dodq2, q2min, q2max, eps);
	return sigmatot;
}

double DM_Particle::Sigma_Electron_Total_Base(double vDM) const
{
	//Numerically integrate the differential cross section
	double q2min						= 0;
	double q2max						= 4.0 * pow(libphysica::Reduced_Mass(mass, mElectron) * vDM, 2.0);
	std::function<double(double)> dodq2 = [this, vDM](double q2) {
		return dSigma_dq2_Electron(sqrt(q2), vDM);
	};
	double eps		= libphysica::Find_Epsilon(dodq2, q2min, q2max, 1.0e-6);
	double sigmatot = libphysica::Integrate(dodq2, q2min, q2max, eps);
	return sigmatot;
}

void DM_Particle::Print_Summary_Base(int MPI_rank) const
{
	if(MPI_rank == 0)
	{
		std::cout << std::endl
				  << "----------------------------------------" << std::endl
				  << "DM particle summary:" << std::endl;

		double massunit			= (mass < keV) ? eV : ((mass < MeV) ? keV : ((mass < GeV) ? MeV : GeV));
		std::string massunitstr = (mass < keV) ? "eV" : ((mass < MeV) ? "keV" : ((mass < GeV) ? "MeV" : "GeV"));
		std::cout << "\tMass:\t\t\t" << In_Units(mass, massunit) << " " << massunitstr << std::endl
				  << "\tSpin:\t\t\t" << spin << std::endl
				  << "\tLow mass:\t\t" << ((low_mass) ? "[x]" : "[ ]") << std::endl;
	}
}

double DM_Particle::dSigma_dER_Nucleus(double ER, const Isotope& target, double vDM) const
{
	double q = sqrt(2.0 * target.mass * ER);
	return 2.0 * target.mass * dSigma_dq2_Nucleus(q, target, vDM);
}

double DM_Particle::Sigma_Total_Nucleus(const Isotope& target, double vDM) const
{
	return Sigma_Nucleus_Total_Base(target, vDM);
}
double DM_Particle::Sigma_Total_Electron(double vDM) const
{
	return Sigma_Electron_Total_Base(vDM);
}

void DM_Particle::Print_Summary(int MPI_rank) const
{
	Print_Summary_Base(MPI_rank);
}

// Scattering angle functions
double DM_Particle::PDF_Scattering_Angle_Nucleus_Base(double cos_alpha, const Isotope& target, double vDM)
{
	double q		= libphysica::Reduced_Mass(target.mass, mass) * vDM * sqrt(2.0 * (1.0 - cos_alpha));
	double q2max	= 4.0 * libphysica::Reduced_Mass(target.mass, mass) * libphysica::Reduced_Mass(target.mass, mass) * vDM * vDM;
	double SigmaTot = Sigma_Total_Nucleus(target, vDM);
	if(SigmaTot != 0.0)
		return q2max / 2.0 / SigmaTot * dSigma_dq2_Nucleus(q, target, vDM);
	else
		return 0.0;
}

double DM_Particle::PDF_Scattering_Angle_Electron_Base(double cos_alpha, double vDM)
{
	double q		= libphysica::Reduced_Mass(mElectron, mass) * vDM * sqrt(2.0 * (1.0 - cos_alpha));
	double q2max	= 4.0 * libphysica::Reduced_Mass(mElectron, mass) * libphysica::Reduced_Mass(mElectron, mass) * vDM * vDM;
	double SigmaTot = Sigma_Total_Electron(vDM);
	if(SigmaTot != 0.0)
		return q2max / 2.0 / SigmaTot * dSigma_dq2_Electron(q, vDM);
	else
		return 0.0;
}

double DM_Particle::CDF_Scattering_Angle_Nucleus_Base(double cos_alpha, const Isotope& target, double vDM)
{
	auto integrand = std::bind(&DM_Particle::PDF_Scattering_Angle_Nucleus_Base, this, std::placeholders::_1, target, vDM);
	double epsilon = libphysica::Find_Epsilon(integrand, -1.0, cos_alpha, 1e-6);
	double cdf	   = libphysica::Integrate(integrand, -1.0, cos_alpha, epsilon);
	return cdf;
}

double DM_Particle::CDF_Scattering_Angle_Electron_Base(double cos_alpha, double vDM)
{
	auto integrand = std::bind(&DM_Particle::PDF_Scattering_Angle_Electron_Base, this, std::placeholders::_1, vDM);
	double epsilon = libphysica::Find_Epsilon(integrand, -1.0, cos_alpha, 1e-6);
	double cdf	   = libphysica::Integrate(integrand, -1.0, cos_alpha, epsilon);
	return cdf;
}

double DM_Particle::Sample_Scattering_Angle_Nucleus_Base(const Isotope& target, double vDM, std::mt19937& PRNG)
{
	double xi						  = libphysica::Sample_Uniform(PRNG, 0.0, 1.0);
	std::function<double(double)> cdf = [this, xi, &target, vDM](double cosa) {
		return xi - CDF_Scattering_Angle_Nucleus(cosa, target, vDM);
	};
	double cos_alpha = libphysica::Find_Root(cdf, -1.0, 1.0, 1e-6);
	return cos_alpha;
}

double DM_Particle::Sample_Scattering_Angle_Electron_Base(double vDM, std::mt19937& PRNG)
{
	double xi						  = libphysica::Sample_Uniform(PRNG, 0.0, 1.0);
	std::function<double(double)> cdf = [this, xi, vDM](double cosa) {
		return xi - CDF_Scattering_Angle_Electron(cosa, vDM);
	};
	double cos_alpha = libphysica::Find_Root(cdf, -1.0, 1.0, 1e-6);
	return cos_alpha;
}

double DM_Particle::PDF_Scattering_Angle_Nucleus(double cos_alpha, const Isotope& target, double vDM)
{
	return PDF_Scattering_Angle_Nucleus_Base(cos_alpha, target, vDM);
}

double DM_Particle::PDF_Scattering_Angle_Electron(double cos_alpha, double vDM)
{
	return PDF_Scattering_Angle_Electron_Base(cos_alpha, vDM);
}

double DM_Particle::CDF_Scattering_Angle_Nucleus(double cos_alpha, const Isotope& target, double vDM)
{
	return CDF_Scattering_Angle_Nucleus_Base(cos_alpha, target, vDM);
}

double DM_Particle::CDF_Scattering_Angle_Electron(double cos_alpha, double vDM)
{
	return CDF_Scattering_Angle_Electron_Base(cos_alpha, vDM);
}

double DM_Particle::Sample_Scattering_Angle_Nucleus(const Isotope& target, double vDM, std::mt19937& PRNG)
{
	return Sample_Scattering_Angle_Nucleus_Base(target, vDM, PRNG);
}

double DM_Particle::Sample_Scattering_Angle_Electron(double vDM, std::mt19937& PRNG)
{
	return Sample_Scattering_Angle_Electron_Base(vDM, PRNG);
}

}	// namespace obscura