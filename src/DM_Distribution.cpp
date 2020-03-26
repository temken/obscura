#include "DM_Distribution.hpp"

#include <iostream>
#include <functional>

#include "Numerics.hpp"
#include "Natural_Units.hpp"

//1. Abstract base class for DM distributions that can be used to compute direct detection recoil spectra.

	//Constructors:
	DM_Distribution::DM_Distribution()
	: name("DM base distribution"), DM_density(0.0), v_domain(std::vector<double>{0.0, 1.0})
	{
	}
	DM_Distribution::DM_Distribution(std::string label, double rhoDM, double vMax, double vMin)
	: name(label), DM_density(rhoDM), v_domain(std::vector<double>{vMin, vMax})
	{
	}

	double DM_Distribution::CDF_Speed(double v) const
	{
		if(v < v_domain[0]) return 0.0;
		else if(v > v_domain[1]) return 1.0;
		else
		{
			std::function<double(double)> integrand = [this] (double v)
			{
				return PDF_Speed(v);
			};
			double cdf = Integrate(integrand, v_domain[0], v, 1.0e-6);
			return cdf;
		}
	}

	Vector DM_Distribution::Average_Velocity() const
	{
		Vector v_average(3);
		// Todo
		return v_average;
	}

	double DM_Distribution::Average_Speed() const
	{
		std::function<double(double)> integrand = [this] (double v)
		{
			return v * PDF_Speed(v);
		};
		double v_average = Integrate(integrand, v_domain[0], v_domain[1], 1.0e-3*km/sec);
		return v_average;
	}

	double DM_Distribution::Eta_Function(double vMin) const
	{
		if(vMin < v_domain[0])
		{
			std::cerr<<"Error in DM_Distribution::Eta_Function(double): vMin = "<<In_Units(vMin,km/sec)<<"km/sec lies below the domain ["<<In_Units(v_domain[0],km/sec)<<"km/sec,"<<In_Units(v_domain[1],km/sec)<<"km/sec]."<<std::endl;
			std::exit(EXIT_FAILURE);
		}
		else if (vMin > v_domain[1])
		{
			return 0.0;
		}
		else
		{
			std::function<double(double)> integrand = [this] (double v)
			{
				return 1.0 / v * PDF_Speed(v);
			};
			double eps = Integrate(integrand, vMin, v_domain[1], 1.0e-5);
			double eta = Integrate(integrand, vMin, v_domain[1], eps);
			return eta;
		}
		
	}

	void DM_Distribution::Print_Summary_Base() const
	{
		std::cout<<"Dark matter distribution - Summary" <<std::endl
				<<"\t"<<name <<std::endl<<std::endl
				<<"\tLocal DM density[GeV/cm^3]:\t" <<In_Units(DM_density,GeV/cm/cm/cm)<<std::endl
				<<"\tSpeed domain [km/sec]:\t\t[" <<In_Units(v_domain[0],km/sec)<<","<<In_Units(v_domain[1],km/sec)<<"]"<<std::endl;
	}