#include "obscura/DM_Halo_Models.hpp"

#include <cmath>

#include "libphysica/Integration.hpp"
#include "libphysica/Natural_Units.hpp"
#include "libphysica/Special_Functions.hpp"
#include "libphysica/Utilities.hpp"

#include "obscura/Astronomy.hpp"

namespace obscura
{
using namespace libphysica::natural_units;

//2. Standard halo model (SHM)
//Constructors:
Standard_Halo_Model::Standard_Halo_Model()
: DM_Distribution("Standard halo model (SHM)", 0.4 * GeV / cm / cm / cm, 0.0, (544.0 + 232.58) * km / sec), v_0(220.0 * km / sec), v_esc(544.0 * km / sec)
{
	vel_observer = libphysica::Vector({0, 220.0 * km / sec, 0}) + libphysica::Vector({11.1 * km / sec, 12.2 * km / sec, 7.3 * km / sec});
	v_observer	 = vel_observer.Norm();
	Normalize_PDF();
	DD_use_eta_function = true;
}

Standard_Halo_Model::Standard_Halo_Model(double rho, double v0, double vobs, double vesc)
: DM_Distribution("Standard halo model (SHM)", rho, 0.0, vesc + vobs), v_0(v0), v_esc(vesc), v_observer(vobs)
{
	vel_observer = libphysica::Vector({0, vobs, 0});
	Normalize_PDF();
	DD_use_eta_function = true;
}

Standard_Halo_Model::Standard_Halo_Model(double rho, double v0, libphysica::Vector& vel_obs, double vesc)
: DM_Distribution("Standard halo model", rho, 0.0, vesc + vel_obs.Norm()), v_0(v0), v_esc(vesc), vel_observer(vel_obs)
{
	v_observer = vel_observer.Norm();
	Normalize_PDF();
	DD_use_eta_function = true;
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
void Standard_Halo_Model::Set_Observer_Velocity(const libphysica::Vector& vel_obs)
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

libphysica::Vector Standard_Halo_Model::Get_Observer_Velocity() const
{
	return vel_observer;
}

//Compute N_esc
void Standard_Halo_Model::Normalize_PDF()
{
	N_esc = erf(v_esc / v_0) - 2 * v_esc / v_0 / sqrt(M_PI) * exp(-v_esc * v_esc / v_0 / v_0);
}

//Distribution functions
double Standard_Halo_Model::PDF_Velocity_SHM(libphysica::Vector vel)
{
	double v = vel.Norm();
	if(v > v_esc || v < v_domain[0])
		return 0.0;
	else
		return 1.0 / N_esc * pow(v_0 * sqrt(M_PI), -3.0) * exp(-1.0 * vel * vel / v_0 / v_0);
}

double Standard_Halo_Model::PDF_Speed_SHM(double v)
{
	if(v < v_domain[0] || v > v_domain[1])
		return 0.0;
	else if(v_observer > 0)
		return v / N_esc / v_0 / sqrt(M_PI) / v_observer * (2 * exp(-(v * v + v_observer * v_observer) / v_0 / v_0) * sinh(2 * v * v_observer / v_0 / v_0) + (exp(-pow(v + v_observer, 2.0) / v_0 / v_0) - exp(-v_esc * v_esc / v_0 / v_0)) * libphysica::StepFunction(fabs(v + v_observer) - v_esc));	 // - (exp(-pow(v - v_observer, 2.0) / v_0 / v_0) - exp(-v_esc * v_esc / v_0 / v_0)) * libphysica::StepFunction(fabs(v - v_observer) - v_esc));
	else
		return 4.0 * v * v / N_esc / sqrt(M_PI) / v_0 / v_0 / v_0 * exp(-v * v / v_0 / v_0) * libphysica::StepFunction(Maximum_DM_Speed() - v);
}

double Standard_Halo_Model::CDF_Speed_SHM(double v)
{
	if(v <= v_domain[0])
		return 0.0;
	else if(v >= v_domain[1])
		return 1.0;
	else if(v_observer > 0)
	{
		double cdf = 1.0 / 2.0 / N_esc * (erf((v - v_observer) / v_0) + erf((v + v_observer) / v_0) - v_0 / v_observer / sqrt(M_PI) * (exp(-(v - v_observer) * (v - v_observer) / v_0 / v_0) - exp(-(v + v_observer) * (v + v_observer) / v_0 / v_0)));
		if(v > (v_esc - v_observer))
			cdf += 1.0 / 2.0 / N_esc * (erf(v_esc / v_0) - erf((v + v_observer) / v_0) - v_0 / v_observer / sqrt(M_PI) * exp(-(v + v_observer) * (v + v_observer) / v_0 / v_0) + (v_0 * v_0 - v * v + (v_esc - v_observer) * (v_esc - v_observer)) / v_observer / v_0 / sqrt(M_PI) * exp(-v_esc * v_esc / v_0 / v_0));
		return cdf;
	}
	else
		return 1.0 / N_esc * (erf(v / v_0) - 2.0 * v / sqrt(M_PI) / v_0 * exp(-v * v / v_0 / v_0));
}

double Standard_Halo_Model::PDF_Velocity(libphysica::Vector vel)
{
	return PDF_Velocity_SHM(vel + vel_observer);
}
double Standard_Halo_Model::PDF_Speed(double v)
{
	return PDF_Speed_SHM(v);
}
double Standard_Halo_Model::CDF_Speed(double v)
{
	return CDF_Speed_SHM(v);
}

//Eta-function for direct detection
double Standard_Halo_Model::Eta_Function_SHM(double vMin)
{
	double xMin = vMin / v_0;
	double xEsc = v_esc / v_0;
	double xE	= v_observer / v_0;
	if(xMin > (xE + xEsc))
		return 0.0;
	else if(fabs(xMin - xE - xEsc) < 1e-8)
		return 0.0;
	else if(xE < 1e-8)
		return 2.0 / N_esc / sqrt(M_PI) / v_0 * (exp(-xMin * xMin) - exp(-xEsc * xEsc));
	else if(xMin > fabs(xE - xEsc))
		return 1.0 / v_0 / 2.0 / N_esc / xE * (erf(xEsc) - erf(xMin - xE) - 2.0 / sqrt(M_PI) * (xE + xEsc - xMin) * exp(-xEsc * xEsc));
	else if(xEsc > xE)
		return 1.0 / v_0 / 2.0 / N_esc / xE * (erf(xMin + xE) - erf(xMin - xE) - 4.0 / sqrt(M_PI) * xE * exp(-xEsc * xEsc));
	else
		return 1.0 / v_0 / xE;
}

double Standard_Halo_Model::Eta_Function(double vMin)
{
	return Eta_Function_SHM(vMin);
}

void Standard_Halo_Model::Print_Summary_SHM()
{
	std::cout << "\tSpeed dispersion v_0[km/sec]:\t" << In_Units(v_0, km / sec) << std::endl
			  << "\tGal. escape velocity [km/sec]:\t" << In_Units(v_esc, km / sec) << std::endl
			  << "\tObserver's velocity [km/sec]:\t" << libphysica::Round(In_Units(vel_observer, km / sec)) << std::endl
			  << "\tObserver's speed [km/sec]:\t" << libphysica::Round(In_Units(v_observer, km / sec)) << std::endl
			  << std::endl;
}

void Standard_Halo_Model::Print_Summary(int mpi_rank)
{
	if(mpi_rank == 0)
	{
		Print_Summary_Base();
		Print_Summary_SHM();
	}
}

//3. Standard halo model++ (SHM++) as proposed by Evans, O'Hare and McCabe [arXiv:1810.11468]
void SHM_Plus_Plus::Compute_Sigmas(double beta)
{
	sigma_r		= sqrt(3.0 * v_0 * v_0 / 2.0 / (3.0 - 2.0 * beta));
	sigma_theta = sqrt(3.0 * v_0 * v_0 * (1.0 - beta) / 2.0 / (3.0 - 2.0 * beta));
	sigma_phi	= sigma_theta;
}

void SHM_Plus_Plus::Normalize_PDF()
{
	N_esc	= erf(v_esc / v_0) - 2 * v_esc / v_0 / sqrt(M_PI) * exp(-v_esc * v_esc / v_0 / v_0);
	N_esc_S = erf(v_esc / sqrt(2.0) / sigma_r) - sqrt((1.0 - beta) / beta) * exp(-v_esc * v_esc / 2.0 / sigma_theta / sigma_theta) * libphysica::Erfi(v_esc / sqrt(2.0) / sigma_r * sqrt(beta / (1.0 - beta)));
}

double SHM_Plus_Plus::PDF_Velocity_S(libphysica::Vector vel)
{
	double v = vel.Norm();
	if(v > v_esc)
		return 0.0;
	else
	{
		double v_r	   = vel[0];
		double v_theta = vel[2];
		double v_phi   = vel[1];
		return 1.0 / N_esc_S / std::pow(2.0 * M_PI, 1.5) / sigma_r / sigma_theta / sigma_phi * exp(-v_r * v_r / 2.0 / sigma_r / sigma_r - v_theta * v_theta / 2.0 / sigma_theta / sigma_theta - v_phi * v_phi / 2.0 / sigma_phi / sigma_phi);
	}
}

double SHM_Plus_Plus::PDF_Speed_S(double v)
{
	auto integrand = [this, v](double cos_theta, double phi) {
		libphysica::Vector vel = libphysica::Spherical_Coordinates(v, acos(cos_theta), phi);
		return v * v * PDF_Velocity_S(vel + vel_observer);
	};
	return libphysica::Integrate_2D(integrand, -1.0, 1.0, 0.0, 2.0 * M_PI);
}

double SHM_Plus_Plus::CDF_Speed_S(double v)
{
	if(v < v_domain[0])
		return 0.0;
	else if(v > v_domain[1])
		return 1.0;
	else
	{
		std::function<double(double)> integrand = [this](double v) {
			return PDF_Speed_S(v);
		};
		return libphysica::Integrate(integrand, v_domain[0], v);
	}
}

double SHM_Plus_Plus::Eta_Function_S(double vMin)
{
	if(vMin < v_domain[0])
		return 0.0;
	else if(vMin > v_domain[1])
		return 1.0;
	else
	{
		std::function<double(double)> integrand = [this](double v) {
			return PDF_Speed_S(v) / v;
		};
		return libphysica::Integrate(integrand, vMin, v_domain[1]);
	}
}

void SHM_Plus_Plus::Interpolate_Eta_Function_S(int v_points)
{
	std::vector<double> v_list	 = libphysica::Linear_Space(v_domain[0], v_domain[1], v_points);
	std::vector<double> eta_list = {};
	for(auto& v : v_list)
		eta_list.push_back(Eta_Function_S(v));
	eta_interpolation_s = libphysica::Interpolation(v_list, eta_list);
}

void SHM_Plus_Plus::Print_Summary_SHMpp()
{
	std::cout << "\tEta parameter:\t" << eta << std::endl
			  << "\tBeta parameter:\t" << beta << std::endl
			  << std::endl;
}

SHM_Plus_Plus::SHM_Plus_Plus(double e, double b)
: Standard_Halo_Model(0.55 * GeV / cm / cm / cm, 233.0 * km / sec, 240.0 * km / sec, 528.0 * km / sec), eta(e), beta(b)
{
	name = "Refined Standard Halo Model (SHM++)";
	Compute_Sigmas(beta);
	Normalize_PDF();
	Interpolate_Eta_Function_S();
}

SHM_Plus_Plus::SHM_Plus_Plus(double rho, double v0, double vobs, double vesc, double e, double b)
: Standard_Halo_Model(rho, v0, vobs, vesc), eta(e), beta(b)
{
	name = "Refined Standard Halo Model (SHM++)";
	Compute_Sigmas(beta);
	Normalize_PDF();
	Interpolate_Eta_Function_S();
}

SHM_Plus_Plus::SHM_Plus_Plus(double rho, double v0, libphysica::Vector& vel_obs, double vesc, double e, double b)
: Standard_Halo_Model(rho, v0, vel_obs, vesc), eta(e), beta(b)
{
	name = "Refined Standard Halo Model (SHM++)";
	Compute_Sigmas(beta);
	Normalize_PDF();
	Interpolate_Eta_Function_S();
}

void SHM_Plus_Plus::Set_Speed_Dispersion(double v0)
{
	v_0 = v0;
	Compute_Sigmas(beta);
	Normalize_PDF();
}

void SHM_Plus_Plus::Set_Eta(double e)
{
	eta = e;
}

void SHM_Plus_Plus::Set_Beta(double b)
{
	beta = b;
	Compute_Sigmas(beta);
	Normalize_PDF();
}

double SHM_Plus_Plus::PDF_Velocity(libphysica::Vector vel)
{
	return (1.0 - eta) * PDF_Velocity_SHM(vel + vel_observer) + eta * PDF_Velocity_S(vel + vel_observer);
}

double SHM_Plus_Plus::PDF_Speed(double v)
{
	return (1.0 - eta) * PDF_Speed_SHM(v) + eta * PDF_Speed_S(v);
}

double SHM_Plus_Plus::CDF_Speed(double v)
{
	return (1.0 - eta) * CDF_Speed_SHM(v) + eta * CDF_Speed_S(v);
}

double SHM_Plus_Plus::Eta_Function(double vMin)
{
	if(vMin < v_domain[0])
	{
		std::cerr << "Error in obscura::SHM_Plus_Plus::Eta_Function: vMin = " << In_Units(vMin, km / sec) << "km/sec lies below the domain [" << In_Units(v_domain[0], km / sec) << "km/sec," << In_Units(v_domain[1], km / sec) << "km/sec]." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	else if(vMin > v_domain[1])
		return 0.0;
	else
		return (1.0 - eta) * Eta_Function_SHM(vMin) + eta * eta_interpolation_s(vMin);
}

void SHM_Plus_Plus::Print_Summary(int mpi_rank)
{
	if(mpi_rank == 0)
	{
		Print_Summary_Base();
		Print_Summary_SHM();
		Print_Summary_SHMpp();
	}
}

}	// namespace obscura