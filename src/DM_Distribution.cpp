#include "obscura/DM_Distribution.hpp"

#include <functional>
#include <iostream>

#include "libphysica/Integration.hpp"
#include "libphysica/Natural_Units.hpp"
#include "libphysica/Special_Functions.hpp"
#include "libphysica/Utilities.hpp"

namespace obscura
{
using namespace libphysica::natural_units;

// 1. Abstract base class for DM distributions that can be used to compute direct detection recoil spectra.
// Constructors:
DM_Distribution::DM_Distribution()
: name("DM base distribution"), v_domain(std::vector<double> {0.0, 1.0}), DM_density(0.0), DD_use_eta_function(false)
{
}
DM_Distribution::DM_Distribution(std::string label, double rhoDM, double vMin, double vMax)
: name(label), v_domain(std::vector<double> {vMin, vMax}), DM_density(rhoDM), DD_use_eta_function(false)
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

double DM_Distribution::PDF_Speed(double v)
{
	auto integrand = [this, v](double cos_theta, double phi) {
		libphysica::Vector vel = libphysica::Spherical_Coordinates(v, acos(cos_theta), phi);
		return v * v * PDF_Velocity(vel);
	};
	return libphysica::Integrate_2D(integrand, -1.0, 1.0, 0.0, 2.0 * M_PI);
}

double DM_Distribution::CDF_Speed(double v)
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
		return libphysica::Integrate(integrand, v_domain[0], v);
	}
}

double DM_Distribution::PDF_Norm()
{
	auto integrand = [this](double v) {
		return PDF_Speed(v);
	};
	return libphysica::Integrate(integrand, v_domain[0], v_domain[1], "Trapezoidal");
}

double DM_Distribution::Differential_DM_Flux(double v, double mDM)
{
	return DM_density / mDM * v * PDF_Speed(v);
}

double DM_Distribution::Total_DM_Flux(double mDM)
{
	auto dFdv = [this, mDM](double v) {
		return Differential_DM_Flux(v, mDM);
	};
	return libphysica::Integrate(dFdv, v_domain[0], v_domain[1]);
}

libphysica::Vector DM_Distribution::Average_Velocity()
{
	libphysica::Vector v_average(3);
	for(unsigned int i = 0; i < v_average.Size(); i++)
	{
		auto integrand = [this, i](libphysica::Vector vel) {
			return vel[i] * PDF_Velocity(vel);
		};
		v_average[i] = libphysica::Integrate_3D(integrand, v_domain[0], v_domain[1], -1.0, 1.0, 0.0, 2.0 * M_PI);
	}
	return v_average;
}

double DM_Distribution::Average_Speed(double vMin)
{
	// 1. Check the domain.
	bool agerage_over_subdomain = true;
	if(vMin < 0)
	{
		vMin				   = v_domain[0];
		agerage_over_subdomain = false;
	}
	else if(vMin < v_domain[0] || vMin >= v_domain[1])
	{
		std::cerr << "Error in DM_Distribution::Average_Speed(): vMin = " << In_Units(vMin, km / sec) << "lies outside the domain [" << In_Units(v_domain[0], km / sec) << "," << In_Units(v_domain[1], km / sec) << "]" << std::endl;
		std::exit(EXIT_FAILURE);
	}

	// 2. Integrate v*f(v)
	std::function<double(double)> integrand = [this](double v) {
		return v * PDF_Speed(v);
	};
	double v_average = libphysica::Integrate(integrand, vMin, v_domain[1]);

	// 3. Re-normalize in case of a sub-domain
	if(agerage_over_subdomain)
		v_average /= (1.0 - CDF_Speed(vMin));

	return v_average;
}

double DM_Distribution::Eta_Function_Base(double vMin)
{
	if(vMin < v_domain[0])
	{
		std::cerr << "Error in obscura::DM_Distribution::Eta_Function_Base(): vMin = " << In_Units(vMin, km / sec) << "km/sec lies below the domain [" << In_Units(v_domain[0], km / sec) << "km/sec," << In_Units(v_domain[1], km / sec) << "km/sec]." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	else if(vMin >= v_domain[1])
	{
		return 0.0;
	}
	else
	{
		std::function<double(double)> integrand = [this](double v) {
			return 1.0 / v * PDF_Speed(v);
		};
		return libphysica::Integrate(integrand, vMin, v_domain[1]);
	}
}

double DM_Distribution::Eta_Function(double vMin)
{
	return Eta_Function_Base(vMin);
}

void DM_Distribution::Print_Summary_Base()
{
	std::cout << "Dark matter distribution - Summary" << std::endl
			  << "\t" << name << std::endl
			  << std::endl
			  << "\tLocal DM density[GeV/cm^3]:\t" << In_Units(DM_density, GeV / cm / cm / cm) << std::endl
			  << "\tSpeed domain [km/sec]:\t\t[" << libphysica::Round(In_Units(v_domain[0], km / sec)) << "," << libphysica::Round(In_Units(v_domain[1], km / sec)) << "]" << std::endl
			  << "\tAverage DM velocity [km/sec]:\t" << libphysica::Round(In_Units(Average_Velocity(), km / sec)) << std::endl
			  << "\tAverage DM speed [km/sec]:\t" << libphysica::Round(In_Units(Average_Speed(), km / sec)) << std::endl
			  << std::endl;
}

void DM_Distribution::Print_Summary(int mpi_rank)
{
	if(mpi_rank == 0)
		Print_Summary_Base();
}

void DM_Distribution::Export_PDF_Speed(std::string file_path, int v_points, bool log_scale)
{
	auto v_list = log_scale ? libphysica::Log_Space(v_domain[0], v_domain[1], v_points) : libphysica::Linear_Space(v_domain[0], v_domain[1], v_points);
	auto pdf	= [this](double v) {
		   return PDF_Speed(v);
	};
	libphysica::Export_Function(file_path, pdf, v_list, {km / sec, sec / km});
}

void DM_Distribution::Export_Eta_Function(std::string file_path, int v_points, bool log_scale)
{
	auto v_list = log_scale ? libphysica::Log_Space(v_domain[0], v_domain[1], v_points) : libphysica::Linear_Space(v_domain[0], v_domain[1], v_points);
	auto eta	= [this](double v) {
		   return Eta_Function(v);
	};
	libphysica::Export_Function(file_path, eta, v_list, {km / sec, sec / km});
}

// 2. Import a tabulated DM distribution from a file (format v[km/sec] :: f(v) [sec/km])
void Imported_DM_Distribution::Check_Normalization()
{
	double norm = pdf_speed.Integrate(v_domain[0], v_domain[1]);
	if(libphysica::Relative_Difference(norm, 1.0) > 1.0e-3)
		std::cout << "Warning in obscura::Imported_DM_Distribution::Check_Normalization(): Imported pdf is not normalized (norm = " << norm << ")." << std::endl
				  << std::endl;
}

void Imported_DM_Distribution::Interpolate_Eta()
{
	auto v_list = libphysica::Linear_Space(v_domain[0], v_domain[1], 500);
	std::vector<std::vector<double>> eta_table;
	for(auto& v : v_list)
		eta_table.push_back({v, Eta_Function_Base(v)});
	eta_function = libphysica::Interpolation(eta_table);
}

Imported_DM_Distribution::Imported_DM_Distribution(double rho, const std::string& filepath)
: DM_Distribution("Imported DM distribution", rho, 0.0, 1.0), file_path(filepath)
{
	DD_use_eta_function = true;
	auto pdf_table		= libphysica::Import_Table(file_path, {km / sec, sec / km});
	pdf_speed			= libphysica::Interpolation(pdf_table);
	v_domain			= pdf_speed.domain;
	Check_Normalization();
	Interpolate_Eta();
}

Imported_DM_Distribution::Imported_DM_Distribution(std::vector<std::vector<double>>& pdf_table, double rho)
: DM_Distribution("Tabulated DM distribution", rho, 0.0, 1.0), file_path("-")
{
	pdf_speed = libphysica::Interpolation(pdf_table);
	v_domain  = pdf_speed.domain;
	Check_Normalization();
	Interpolate_Eta();
}

double Imported_DM_Distribution::PDF_Speed(double v)
{
	if(v < v_domain[0] || v > v_domain[1])
		return 0.0;
	else
		return pdf_speed(v);
}

double Imported_DM_Distribution::Eta_Function(double vMin)
{
	if(vMin < v_domain[0])
	{
		std::cerr << "Error in obscura::Imported_DM_Distribution::Eta_Function(): vMin = " << In_Units(vMin, km / sec) << "km/sec lies below the domain [" << In_Units(v_domain[0], km / sec) << "km/sec," << In_Units(v_domain[1], km / sec) << "km/sec]." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	else if(vMin > v_domain[1])
		return 0.0;
	else
		return eta_function(vMin);
}

void Imported_DM_Distribution::Print_Summary(int mpi_rank)
{
	if(mpi_rank == 0)
	{
		Print_Summary_Base();
		std::cout << "\tFile path:\t" << file_path << std::endl
				  << std::endl;
	}
}

}	// namespace obscura
