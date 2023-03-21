#include "obscura/Direct_Detection_Crystal.hpp"

#include <cmath>
#include <numeric>

#include "libphysica/Integration.hpp"
#include "libphysica/Natural_Units.hpp"
#include "libphysica/Special_Functions.hpp"
#include "libphysica/Utilities.hpp"

#include "obscura/Target_Atom.hpp"

namespace obscura
{
using namespace libphysica::natural_units;

// 1. Event spectra and rates

double Minimum_Electron_Energy(int Q, const Crystal& target)
{
	return target.epsilon * (Q - 1) + target.energy_gap;
}

unsigned int Electron_Hole_Pairs(double Ee, const Crystal& target)
{
	return std::floor((Ee - target.energy_gap) / target.epsilon + 1);
}

bool dRdE_Crystal_warned = false;
double dRdEe_Crystal(double Ee, const DM_Particle& DM, DM_Distribution& DM_distr, Crystal& target_crystal)
{
	if(Ee > target_crystal.E_max)
	{
		if(!dRdE_Crystal_warned)
		{
			std::cerr << libphysica::Formatted_String("Warning", "Yellow", true) << " in dRdEe_Crystal: Ee lies beyond the tabulated crystal form factor. Return 0." << std::endl
					  << "\tEe = " << libphysica::Round(Ee / eV) << " eV > E_max = " << libphysica::Round(target_crystal.E_max / eV) << " eV" << std::endl
					  << "\t(Warning will not be repeated.)" << std::endl;
			dRdE_Crystal_warned = true;
		}
		return 0;
	}
	double N_T		= 1.0 / target_crystal.M_cell;
	double integral = 0.0;
	for(int qi = 0; qi < target_crystal.N_q; qi++)
	{
		double q	= (qi + 1) * target_crystal.dq;
		double vMin = vMinimal_Electrons(q, Ee, DM.mass);
		double vMax = DM_distr.Maximum_DM_Speed();
		if(vMin > vMax)
			continue;
		else if(DM.DD_use_eta_function && DM_distr.DD_use_eta_function)
		{
			double vDM = 1e-3;	 // cancels in v^2 * dSigma/dq^2
			integral += 2.0 * q * target_crystal.dq * DM_distr.DM_density / DM.mass * DM_distr.Eta_Function(vMin) * vDM * vDM * DM.d2Sigma_dq2_dEe_Crystal(q, Ee, vDM, target_crystal);
		}
		else
		{
			auto integrand = [&DM_distr, &DM, q, Ee, &target_crystal](double v) {
				return DM_distr.Differential_DM_Flux(v, DM.mass) * DM.d2Sigma_dq2_dEe_Crystal(q, Ee, v, target_crystal);
			};
			integral += 2.0 * q * target_crystal.dq * libphysica::Integrate(integrand, vMin, vMax);
		}
	}
	return N_T * integral;
}

double R_Q_Crystal(int Q, const DM_Particle& DM, DM_Distribution& DM_distr, Crystal& target_crystal)
{
	// Energy threshold
	double Emin = Minimum_Electron_Energy(Q, target_crystal);
	double Emax = Minimum_Electron_Energy(Q + 1, target_crystal);
	// Integrate over energies
	double sum = 0.0;
	for(int Ei = (Emin / target_crystal.dE); Ei < target_crystal.N_E; Ei++)
	{
		double E = (Ei + 1) * target_crystal.dE;
		if(E > Emax)
			break;
		sum += target_crystal.dE * dRdEe_Crystal(E, DM, DM_distr, target_crystal);
	}
	return sum;
}

double R_total_Crystal(int Qthreshold, const DM_Particle& DM, DM_Distribution& DM_distr, Crystal& target_crystal)
{
	// Energy threshold
	double E_min = Minimum_Electron_Energy(Qthreshold, target_crystal);
	// Integrate over energies
	double sum = 0.0;
	for(int Ei = (E_min / target_crystal.dE); Ei < target_crystal.N_E; Ei++)
	{
		double E = (Ei + 1) * target_crystal.dE;
		sum += target_crystal.dE * dRdEe_Crystal(E, DM, DM_distr, target_crystal);
	}
	return sum;
}

// 2. Electron recoil direct detection experiment with semiconductor target
DM_Detector_Crystal::DM_Detector_Crystal()
: DM_Detector("Crystal experiment", gram * year, "Electrons"), target_crystal(Crystal("Si"))
{
}

DM_Detector_Crystal::DM_Detector_Crystal(std::string label, double expo, std::string crys)
: DM_Detector(label, expo, "Electrons"), target_crystal(Crystal(crys))
{
}

// DM functions
double DM_Detector_Crystal::Maximum_Energy_Deposit(DM_Particle& DM, const DM_Distribution& DM_distr) const
{
	return DM.mass / 2.0 * pow(DM_distr.Maximum_DM_Speed(), 2.0);
}

double DM_Detector_Crystal::Minimum_DM_Mass(DM_Particle& DM, const DM_Distribution& DM_distr) const
{
	return 2.0 * energy_threshold * pow(DM_distr.Maximum_DM_Speed(), -2.0);
}

double DM_Detector_Crystal::Minimum_DM_Speed(DM_Particle& DM) const
{
	return sqrt(2.0 * energy_threshold / DM.mass);
}

double DM_Detector_Crystal::dRdE(double E, const DM_Particle& DM, DM_Distribution& DM_distr)
{
	return flat_efficiency * dRdEe_Crystal(E, DM, DM_distr, target_crystal);
}

double DM_Detector_Crystal::DM_Signals_Total(const DM_Particle& DM, DM_Distribution& DM_distr)
{
	double N = 0;
	if(statistical_analysis == "Binned Poisson")
	{
		std::vector<double> binned_events = DM_Signals_Binned(DM, DM_distr);
		N								  = std::accumulate(binned_events.begin(), binned_events.end(), 0.0);
	}
	else if(using_energy_threshold || statistical_analysis == "Maximum Gap")
	{
		std::function<double(double)> spectrum = [this, &DM, &DM_distr](double E) {
			return dRdE(E, DM, DM_distr);
		};
		N = exposure * libphysica::Integrate(spectrum, energy_threshold, energy_max);
	}
	else if(using_Q_threshold)
	{
		N = exposure * flat_efficiency * R_total_Crystal(Q_threshold, DM, DM_distr, target_crystal);
	}
	return N;
}

std::vector<double> DM_Detector_Crystal::DM_Signals_Binned(const DM_Particle& DM, DM_Distribution& DM_distr)
{
	if(statistical_analysis != "Binned Poisson")
	{
		std::cerr << libphysica::Formatted_String("Error", "Red", true) << " in obscura::DM_Detector_Crystal::DM_Signals_Binned(const DM_Particle&, DM_Distribution&): The statistical analysis is " << statistical_analysis << ", not 'Binned Poisson'." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	else if(using_energy_bins)
	{
		return DM_Signals_Energy_Bins(DM, DM_distr);
	}
	else if(using_Q_bins)
	{
		return DM_Signals_Q_Bins(DM, DM_distr);
	}
	else
	{
		std::cerr << libphysica::Formatted_String("Error", "Red", true) << " in obscura::DM_Detector_Crystal::DM_Signals_Binned(): Statistical analysis is 'Binned Poisson' but no bins have been defined. This should not happen ever." << std::endl;
		std::exit(EXIT_FAILURE);
	}
}

// Q spectrum
void DM_Detector_Crystal::Use_Q_Threshold(unsigned int Q_thr)
{
	Initialize_Poisson();
	using_Q_threshold = true;
	Q_threshold		  = Q_thr;
	energy_threshold  = Minimum_Electron_Energy(Q_threshold, target_crystal);
	energy_max		  = Minimum_Electron_Energy(target_crystal.Q_max, target_crystal);
}

void DM_Detector_Crystal::Use_Q_Bins(unsigned int Q_thr, unsigned int N_bins)
{
	using_Q_bins	   = true;
	Q_threshold		   = Q_thr;
	unsigned int Q_max = (N_bins == 0) ? target_crystal.Q_max : N_bins + Q_threshold - 1;
	N_bins			   = Q_max - Q_threshold + 1;
	Initialize_Binned_Poisson(N_bins);

	energy_threshold = Minimum_Electron_Energy(Q_threshold, target_crystal);
	energy_max		 = Minimum_Electron_Energy(Q_max, target_crystal);
}

std::vector<double> DM_Detector_Crystal::DM_Signals_Q_Bins(const DM_Particle& DM, DM_Distribution& DM_distr)
{
	if(!using_Q_bins)
	{
		std::cerr << libphysica::Formatted_String("Error", "Red", true) << " in obscura::DM_Detector_Crystal::DM_Signals_Q_Bins(const DM_Particle&,DM_Distribution&): Not using Q bins." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	else
	{
		std::vector<double> signals;
		for(unsigned int Q = Q_threshold; Q < Q_threshold + number_of_bins; Q++)
		{
			signals.push_back(exposure * flat_efficiency * bin_efficiencies[Q - 1] * R_Q_Crystal(Q, DM, DM_distr, target_crystal));
		}
		return signals;
	}
}

void DM_Detector_Crystal::Print_Summary(int MPI_rank) const
{
	Print_Summary_Base();
	std::cout << std::endl
			  << "\tElectron recoil experiment (semiconductor)." << std::endl
			  << "\tTarget:\t\t\t" << target_crystal.name << " semiconductor" << std::endl;
	if(using_Q_threshold || using_Q_bins)
		std::cout << "\teh pair threshold:\t" << Q_threshold << std::endl
				  << "----------------------------------------" << std::endl
				  << std::endl;
}

}	// namespace obscura
