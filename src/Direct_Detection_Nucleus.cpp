#include "obscura/Direct_Detection_Nucleus.hpp"

#include <algorithm>   //for std::min_element, std::max_element, std::sort
#include <cmath>
#include <numeric>	 //for std::accumulate

#include "libphysica/Special_Functions.hpp"
#include "libphysica/Statistics.hpp"
#include "libphysica/Utilities.hpp"

namespace obscura
{
using namespace libphysica::natural_units;

//1. Theoretical nuclear recoil spectrum
double dRdER_Nucleus(double ER, const DM_Particle& DM, DM_Distribution& DM_distr, const Isotope& target_isotope)
{
	double vMin = vMinimal_Nucleus(ER, DM.mass, target_isotope.mass);
	double vMax = DM_distr.Maximum_DM_Speed();
	if(vMin > vMax)
		return 0.0;
	else if(DM.DD_use_eta_function && DM_distr.DD_use_eta_function)
	{
		double rhoDM = DM_distr.DM_density * DM.fractional_density;
		double vDM	 = 1.0e-3;	 //cancels when eta function can be used
		return 1.0 / target_isotope.mass * rhoDM / DM.mass * (vDM * vDM * DM.dSigma_dER_Nucleus(ER, target_isotope, vDM)) * DM_distr.Eta_Function(vMin);
	}
	else
	{
		auto integrand = [ER, &DM, &DM_distr, &target_isotope](double v) {
			return DM_distr.Differential_DM_Flux(v, DM.mass) * DM.dSigma_dER_Nucleus(ER, target_isotope, v);
		};
		double integral = libphysica::Integrate(integrand, vMin, vMax);
		return DM.fractional_density / target_isotope.mass * integral;
	}
}

double dRdER_Nucleus(double ER, const DM_Particle& DM, DM_Distribution& DM_distr, const Nucleus& target_nucleus)
{
	double dRate = 0.0;
	for(unsigned int i = 0; i < target_nucleus.Number_of_Isotopes(); i++)
		dRate += target_nucleus[i].abundance * dRdER_Nucleus(ER, DM, DM_distr, target_nucleus[i]);

	return dRate;
}

//2. Nuclear recoil direct detection experiment
//Constructors
DM_Detector_Nucleus::DM_Detector_Nucleus()
: DM_Detector("Nuclear recoil experiment", kg * day, "Nuclei"), target_nuclei({Get_Nucleus(54)}), relative_mass_fractions({1.0}), energy_resolution(0.0), using_efficiency_tables(false)
{
}

DM_Detector_Nucleus::DM_Detector_Nucleus(std::string label, double expo, std::vector<Nucleus> nuclei, std::vector<double> abund)
: DM_Detector(label, expo, "Nuclei"), target_nuclei(nuclei), energy_resolution(0.0), using_efficiency_tables(false)
{
	double tot = std::accumulate(abund.begin(), abund.end(), 0.0);
	if(abund.empty() || tot > 1.0)
	{
		//Compute relative abundance of the target nuclei by weight.
		for(unsigned int i = 0; i < target_nuclei.size(); i++)
		{
			double proportion = (abund.empty()) ? 1.0 : abund[i];
			double mi		  = proportion * target_nuclei[i].Average_Nuclear_Mass();
			relative_mass_fractions.push_back(mi);
		}
		double Mtot = std::accumulate(relative_mass_fractions.begin(), relative_mass_fractions.end(), 0.0);
		for(unsigned int i = 0; i < relative_mass_fractions.size(); i++)
			relative_mass_fractions[i] /= Mtot;
	}
	else
		relative_mass_fractions = abund;
}

double DM_Detector_Nucleus::Maximum_Energy_Deposit(const DM_Particle& DM, const DM_Distribution& DM_distr) const
{
	double vDM	= DM_distr.Maximum_DM_Speed();
	double Emax = 0.0;
	for(unsigned int i = 0; i < target_nuclei.size(); i++)
	{
		for(unsigned int j = 0; j < target_nuclei[i].Number_of_Isotopes(); j++)
		{
			double ERmax = Maximum_Nuclear_Recoil_Energy(vDM, DM.mass, target_nuclei[i][j].mass);
			if(ERmax > Emax && DM.Sigma_Total_Nucleus(target_nuclei[i][j], vDM) > 0.0)
				Emax = ERmax;
		}
	}
	return Emax + 6.0 * energy_resolution;
}

double DM_Detector_Nucleus::Minimum_DM_Mass(DM_Particle& DM, const DM_Distribution& DM_distr) const
{
	std::vector<double> aux;
	double vMax = DM_distr.Maximum_DM_Speed();
	for(unsigned int i = 0; i < target_nuclei.size(); i++)
	{
		for(unsigned int j = 0; j < target_nuclei[i].Number_of_Isotopes(); j++)
		{
			double mMin = target_nuclei[i][j].mass / (sqrt(2.0 * target_nuclei[i][j].mass / (energy_threshold - 2.0 * energy_resolution)) * vMax - 1.0);
			if(DM.Sigma_Total_Nucleus(target_nuclei[i][j], vMax) > 0.0)
				aux.push_back(mMin);
		}
	}
	return *std::min_element(aux.begin(), aux.end());
}

void DM_Detector_Nucleus::Set_Resolution(double res)
{
	energy_resolution = res;
}

void DM_Detector_Nucleus::Import_Efficiency(std::string filename, double dim)
{
	using_efficiency_tables							  = true;
	std::vector<std::vector<double>> efficiency_table = libphysica::Import_Table(filename);
	libphysica::Interpolation eff(efficiency_table, dim);
	efficiencies.push_back(eff);
}

void DM_Detector_Nucleus::Import_Efficiency(std::vector<std::string> filenames, double dim)
{
	efficiencies.clear();
	for(unsigned int i = 0; i < filenames.size(); i++)
		Import_Efficiency(filenames[i], dim);
}

double DM_Detector_Nucleus::dRdE(double E, const DM_Particle& DM, DM_Distribution& DM_distr)
{

	double dR = 0.0;
	if(energy_resolution < 1e-6 * eV)
	{
		double eff = 1.0;
		for(unsigned int i = 0; i < target_nuclei.size(); i++)
		{
			if(using_efficiency_tables)
			{
				if(efficiencies.size() == 1)
					eff = efficiencies[0](E);
				else if(efficiencies.size() == target_nuclei.size())
					eff = efficiencies[i](E);
			}
			dR += eff * flat_efficiency * relative_mass_fractions[i] * dRdER_Nucleus(E, DM, DM_distr, target_nuclei[i]);
		}
	}
	else
	{
		//Find minimum and maximum ER contributing to dR/dE(E):
		std::vector<double> aux = {E - 6.0 * energy_resolution, energy_threshold - 3.0 * energy_resolution, 2.0 * energy_resolution};
		double eMin				= *std::max_element(aux.begin(), aux.end());
		double eMax				= E + 6.0 * energy_resolution;

		//Convolute theoretical spectrum with Gaussian
		std::function<double(double)> integrand = [this, E, &DM, &DM_distr](double ER) {
			double dRtheory = 0.0;
			for(unsigned int i = 0; i < target_nuclei.size(); i++)
			{
				double eff = 1.0;
				if(using_efficiency_tables)
				{
					if(efficiencies.size() == 1)
						eff = efficiencies[0](E);
					else if(efficiencies.size() == target_nuclei.size())
						eff = efficiencies[i](E);
				}
				dRtheory += eff * flat_efficiency * relative_mass_fractions[i] * dRdER_Nucleus(ER, DM, DM_distr, target_nuclei[i]);
			}
			return libphysica::PDF_Gauss(E, ER, energy_resolution) * dRtheory;
		};
		dR = libphysica::Integrate(integrand, eMin, eMax);
	}
	return dR;
}

double DM_Detector_Nucleus::Minimum_DM_Speed(const DM_Particle& DM) const
{
	double Emin = energy_threshold - 2.0 * energy_resolution;
	double vcut = 1.0;
	for(unsigned int i = 0; i < target_nuclei.size(); i++)
	{
		for(unsigned int j = 0; j < target_nuclei[i].Number_of_Isotopes(); j++)
		{
			double vmin = vMinimal_Nucleus(Emin, DM.mass, target_nuclei[i][j].mass);
			if(vmin < vcut && DM.Sigma_Total_Nucleus(target_nuclei[i][j], 1.0e-3) > 0.0)
				vcut = vmin;
		}
	}
	return vcut;
}

void DM_Detector_Nucleus::Print_Summary(int MPI_rank) const
{
	if(MPI_rank == 0)
	{
		Print_Summary_Base(MPI_rank);
		std::cout << std::endl
				  << "\tNuclear recoil experiment." << std::endl
				  << "\tNuclear targets:" << std::endl
				  << "\t\tNucl.\tabund." << std::endl;
		for(unsigned int i = 0; i < target_nuclei.size(); i++)
		{
			std::cout << "\t\t" << target_nuclei[i].name << "\t" << libphysica::Round(100.0 * relative_mass_fractions[i]) << "%" << std::endl;
			// target_nuclei[i].Print_Summary();
		}
		std::cout << "\tThreshold [keV]:\t" << In_Units(energy_threshold, keV) << std::endl
				  << "\tER_max [keV]:\t\t" << In_Units(energy_max, keV) << std::endl
				  << "\tER resolution [keV]:\t" << In_Units(energy_resolution, keV) << std::endl
				  << "----------------------------------------" << std::endl
				  << std::endl;
	}
}

}	// namespace obscura
