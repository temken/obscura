#include "DM_Distribution.hpp"

#include <cmath>
#include <functional>
#include <iostream>

//Headers from libphysica library
#include "Natural_Units.hpp"
#include "Numerics.hpp"

#include "Astronomy.hpp"

namespace obscura
{
using namespace libphysica::natural_units;

//1. Abstract base class for DM distributions that can be used to compute direct detection recoil spectra.
//Constructors:
DM_Distribution::DM_Distribution()
: name("DM base distribution"), v_domain(std::vector<double>{0.0, 1.0}), DM_density(0.0)
{
}
DM_Distribution::DM_Distribution(std::string label, double rhoDM, double vMin, double vMax)
: name(label), v_domain(std::vector<double>{vMin, vMax}), DM_density(rhoDM)
{
}

double DM_Distribution::Minimum_DM_Speed() const
{
	return v_domain[0];
}
double DM_Distribution::Maximum_DM_Speed() const
{
	return v_domain[1];
}

double DM_Distribution::CDF_Speed(double v) const
{
	if(v < v_domain[0])
		return 0.0;
	else if(v > v_domain[1])
		return 1.0;
	else
	{
		std::function<double(double)> integrand = [this](double v) {
			return PDF_Speed(v);
		};
		double cdf = libphysica::Integrate(integrand, v_domain[0], v, 1.0e-6);
		return cdf;
	}
}

libphysica::Vector DM_Distribution::Average_Velocity() const
{
	libphysica::Vector v_average(3);
	// Todo
	return v_average;
}

double DM_Distribution::Average_Speed(double vMin) const
{
	// 1. Check the domain.
	bool agerage_over_subdomain = true;
	if(vMin < 0)
	{
		vMin				   = v_domain[0];
		agerage_over_subdomain = false;
	}
	else if(vMin < v_domain[0] || vMin > v_domain[1])
	{
		std::cerr << "Error in DM_Distribution::Average_Speed(): vMin = " << In_Units(vMin, km / sec) << "lies outside the domain [" << In_Units(v_domain[0], km / sec) << "," << In_Units(v_domain[1], km / sec) << "]" << std::endl;
		std::exit(EXIT_FAILURE);
	}

	// 2. Integrate v*f(v)
	std::function<double(double)> integrand = [this](double v) {
		return v * PDF_Speed(v);
	};
	double v_average = libphysica::Integrate(integrand, vMin, v_domain[1], 1.0e-3 * km / sec);
	
	// 3. Re-normalize in case of a sub-domain
	if(agerage_over_subdomain)
		v_average /= (1.0 - CDF_Speed(vMin));
	
	return v_average;
}

double DM_Distribution::Eta_Function(double vMin) const
{
	if(vMin < v_domain[0])
	{
		std::cerr << "Error in obscura::DM_Distribution::Eta_Function(double): vMin = " << In_Units(vMin, km / sec) << "km/sec lies below the domain [" << In_Units(v_domain[0], km / sec) << "km/sec," << In_Units(v_domain[1], km / sec) << "km/sec]." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	else if(vMin > v_domain[1])
	{
		return 0.0;
	}
	else
	{
		std::function<double(double)> integrand = [this](double v) {
			return 1.0 / v * PDF_Speed(v);
		};
		double eps = libphysica::Integrate(integrand, vMin, v_domain[1], 1.0e-5);
		double eta = libphysica::Integrate(integrand, vMin, v_domain[1], eps);
		return eta;
	}
}

void DM_Distribution::Print_Summary_Base(int MPI_rank) const
{
	if(MPI_rank == 0)
	{
		std::cout << "Dark matter distribution - Summary" << std::endl
				  << "\t" << name << std::endl
				  << std::endl
				  << "\tLocal DM density[GeV/cm^3]:\t" << In_Units(DM_density, GeV / cm / cm / cm) << std::endl
				  << "\tSpeed domain [km/sec]:\t\t[" << In_Units(v_domain[0], km / sec) << "," << In_Units(v_domain[1], km / sec) << "]" << std::endl
				  << "\tAverage DM velocity [km/sec]:\t" << In_Units(Average_Velocity(), km / sec) << std::endl
				  << "\tAverage DM speed [km/sec]:\t" << libphysica::Round(In_Units(Average_Speed(), km / sec)) << std::endl
				  << std::endl;
	}
}

//2. Standard halo model (SHM)
//Constructors:
Standard_Halo_Model::Standard_Halo_Model()
: DM_Distribution("Standard halo model (SHM)", 0.4 * GeV / cm / cm / cm, 0.0, (544.0 + 232.58) * km / sec), v_0(220.0 * km / sec), v_esc(544.0 * km / sec)
{
	vel_observer = libphysica::Vector({0, 220.0 * km / sec, 0}) + libphysica::Vector({11.1 * km / sec, 12.2 * km / sec, 7.3 * km / sec});
	v_observer	 = vel_observer.Norm();
	Normalize_PDF();
}

Standard_Halo_Model::Standard_Halo_Model(double rho, double v0, double vobs, double vesc)
: DM_Distribution("Standard halo model (SHM)", rho, 0.0, vesc + vobs), v_0(v0), v_esc(vesc), v_observer(vobs)
{
	vel_observer = libphysica::Vector({0, vobs, 0});
	Normalize_PDF();
}

Standard_Halo_Model::Standard_Halo_Model(double rho, double v0, libphysica::Vector& vel_obs, double vesc)
: DM_Distribution("Standard halo model", rho, 0.0, vesc + vel_obs.Norm()), v_0(v0), v_esc(vesc), vel_observer(vel_obs)
{
	v_observer = vel_observer.Norm();
	Normalize_PDF();
}

//Set SHM parameters
void Standard_Halo_Model::Set_Speed_Dispersion(double v0)
{
	v_0 = v0;
	Normalize_PDF();
}
void Standard_Halo_Model::Set_Escape_Velocity(double vesc)
{
	v_esc = vesc;

	v_domain[1] = vesc + v_observer;
	Normalize_PDF();
}
void Standard_Halo_Model::Set_Observer_Velocity(libphysica::Vector& vel_obs)
{
	vel_observer = vel_obs;
	v_observer	 = vel_observer.Norm();

	v_domain[1] = v_esc + v_observer;
}
void Standard_Halo_Model::Set_Observer_Velocity(int day, int month, int year, int hour, int minute)
{
	double nJ2000 = Fractional_Days_since_J2000(day, month, year, hour, minute);
	vel_observer  = Earth_Velocity(nJ2000);
	v_observer	  = vel_observer.Norm();

	v_domain[1] = v_esc + v_observer;
}

//Compute N_esc
void Standard_Halo_Model::Normalize_PDF()
{
	N_esc = erf(v_esc / v_0) - 2 * v_esc / v_0 / sqrt(M_PI) * exp(-v_esc * v_esc / v_0 / v_0);
}

//Distribution functions
double Standard_Halo_Model::PDF_Velocity(libphysica::Vector vel) const
{
	double v = vel.Norm();
	if(v > v_domain[1] || v < v_domain[0])
		return 0.0;
	else
	{
		return 1.0 / N_esc * pow(v_0 * sqrt(M_PI), -3.0) * exp(-1.0 * (vel + vel_observer) * (vel + vel_observer) / v_0 / v_0);
	}
}
double Standard_Halo_Model::PDF_Speed(double v) const
{
	return v / N_esc / v_0 / sqrt(M_PI) / v_observer * (2 * exp(-(v * v + v_observer * v_observer) / v_0 / v_0) * sinh(2 * v * v_observer / v_0 / v_0) + (exp(-pow(v + v_observer, 2.0) / v_0 / v_0) - exp(-v_esc * v_esc / v_0 / v_0)) * libphysica::StepFunction(abs(v + v_observer) - v_esc) - (exp(-pow(v - v_observer, 2.0) / v_0 / v_0) - exp(-v_esc * v_esc / v_0 / v_0)) * libphysica::StepFunction(abs(v - v_observer) - v_esc));
}

//Eta-function for direct detection
double Standard_Halo_Model::Eta_Function(double vMin) const
{
	double xMin = vMin / v_0;
	double xEsc = v_esc / v_0;
	double xE	= v_observer / v_0;
	if(xMin > (xE + xEsc))
		return 0.0;
	else if(fabs(xMin - xE - xEsc) < 1e-8)
		return 0.0;
	else if(xMin > fabs(xE - xEsc))
		return 1.0 / v_0 / 2.0 / N_esc / xE * (erf(xEsc) - erf(xMin - xE) - 2.0 / sqrt(M_PI) * (xE + xEsc - xMin) * exp(-xEsc * xEsc));
	else if(xEsc > xE)
		return 1.0 / v_0 / 2.0 / N_esc / xE * (erf(xMin + xE) - erf(xMin - xE) - 4.0 / sqrt(M_PI) * xE * exp(-xEsc * xEsc));
	else
		return 1.0 / v_0 / xE;
}

void Standard_Halo_Model::Print_Summary(int MPI_rank) const
{
	if(MPI_rank == 0)
	{
		Print_Summary_Base(MPI_rank);
		std::cout << "\tSpeed dispersion v_0[km/sec]:\t" << In_Units(v_0, km / sec) << std::endl
				  << "\tGal. escape velocity [km/sec]:\t" << In_Units(v_esc, km / sec) << std::endl
				  << "\tObserver's velocity [km/sec]:\t" << In_Units(vel_observer, km / sec) << std::endl
				  << "\tObserver's speed [km/sec]:\t" << In_Units(v_observer, km / sec) << std::endl
				  << std::endl;
	}
}

}	// namespace obscura
