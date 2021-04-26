#include "obscura/Direct_Detection_Migdal.hpp"

#include <cmath>

#include "libphysica/Natural_Units.hpp"
#include "libphysica/Numerics.hpp"
#include "libphysica/Statistics.hpp"
#include "libphysica/Utilities.hpp"

#include "obscura/Direct_Detection_Ionization.hpp"

namespace obscura
{
using namespace libphysica::natural_units;

//1. Event spectra and rates
double dRdEe_Migdal(double Ee, const DM_Particle& DM, DM_Distribution& DM_distr, Isotope& isotope, Atomic_Electron& shell)
{
	std::function<double(double)> ER_integrand = [Ee, &DM, &DM_distr, &isotope, &shell](double ER) {
		double q	= sqrt(2.0 * isotope.mass * ER);
		double qe	= mElectron / isotope.mass * q;
		double vMin = (q > 0) ? vMinimal_Electrons(q, (Ee + shell.binding_energy), DM.mass) : 1.0e-20 * cm / sec;
		// check eq 5 of https://arxiv.org/pdf/1711.09906.pdf
		if(vMin >= DM_distr.Maximum_DM_Speed())
			return 0.0;
		if(DM.DD_use_eta_function && DM_distr.DD_use_eta_function)
		{
			double vDM	  = 1.0e-3;	  // cancels
			double result = DM_distr.DM_density / DM.mass * DM_distr.Eta_Function(vMin) * DM.dSigma_dER_Nucleus(ER, isotope, vDM) * vDM * vDM * shell.Ionization_Form_Factor(qe, Ee);
			if(result == 0.0)
				std::cout << vMin / km * sec << "\t" << DM_distr.Eta_Function(vMin) << "\t" << DM.dSigma_dER_Nucleus(ER, isotope, vDM) * vDM * vDM << "\t" << shell.Ionization_Form_Factor(qe, Ee) << std::endl;
			return result;
		}
		else
		{
			std::function<double(double)> v_integrand = [ER, &DM_distr, &DM, &isotope](double v) {
				return DM_distr.Differential_DM_Flux(v, DM.mass) * DM.dSigma_dER_Nucleus(ER, isotope, v);
			};
			double eps		  = libphysica::Find_Epsilon(v_integrand, vMin, DM_distr.Maximum_DM_Speed(), 1.0e-4);
			double v_integral = libphysica::Integrate(v_integrand, vMin, DM_distr.Maximum_DM_Speed(), eps);
			return shell.Ionization_Form_Factor(qe, Ee) * v_integral;
		}
	};
	double vMax = DM_distr.Maximum_DM_Speed();

	double qMin	  = shell.binding_energy / vMax;
	double qMax	  = 2.0 * libphysica::Reduced_Mass(DM.mass, isotope.mass) * vMax;
	double ER_min = qMin * qMin / 2.0 / isotope.mass;
	double ER_max = qMax * qMax / 2.0 / isotope.mass;

	std::vector<double> ER_list = libphysica::Log_Space(ER_min, ER_max, 50);
	double dln_ER				= log(ER_list[1] / ER_list[0]);
	double ER_integral			= 0.0;
	for(auto& ER : ER_list)
		ER_integral += dln_ER * ER * ER_integrand(ER);
	return 1.0 / isotope.mass / 4.0 / Ee * ER_integral;
}

double dRdEe_Migdal(double Ee, const DM_Particle& DM, DM_Distribution& DM_distr, Element& element, Atomic_Electron& shell)
{
	double result = 0.0;
	for(auto& isotope : element.isotopes)
		result += isotope.abundance * dRdEe_Migdal(Ee, DM, DM_distr, isotope, shell);
	return result;
}

double R_ne_Migdal(unsigned int ne, const DM_Particle& DM, DM_Distribution& DM_distr, Element& element, Atomic_Electron& shell)
{
	double sum = 0.0;
	for(unsigned int ki = 0; ki < shell.Nk; ki++)
	{
		double k  = shell.k_Grid[ki];
		double Ee = k * k / 2.0 / mElectron;
		sum += log(10.0) * shell.dlogk * k * k / mElectron * PDF_ne(ne, Ee, shell) * dRdEe_Migdal(Ee, DM, DM_distr, element, shell);
	}
	return sum;
}

double R_ne_Migdal(unsigned int ne, const DM_Particle& DM, DM_Distribution& DM_distr, Element& element, Atom& atom)
{
	double result = 0.0;
	for(unsigned int i = 0; i < atom.electrons.size(); i++)
		for(auto& electron : atom.electrons)
			result += R_ne_Migdal(ne, DM, DM_distr, element, electron);
	return result;
}

double R_PE_Migdal_aux(unsigned int nPE, double mu_PE, double sigma_PE, std::vector<double> R_ne_spectrum)
{
	double sum = 0.0;
	for(int ne = 1; ne < 16; ne++)
		sum += libphysica::PDF_Gauss(nPE, mu_PE * ne, sqrt(ne) * sigma_PE) * R_ne_spectrum[ne - 1];
	return sum;
}

double R_PE_Migdal(unsigned int nPE, double mu_PE, double sigma_PE, const DM_Particle& DM, DM_Distribution& DM_distr, Element& element, Atomic_Electron& shell, std::vector<double> electron_spectrum)
{
	if(electron_spectrum.empty())
		for(unsigned ne = 1; ne < 16; ne++)
			electron_spectrum.push_back(R_ne_Migdal(ne, DM, DM_distr, element, shell));
	return R_PE_Migdal_aux(nPE, mu_PE, sigma_PE, electron_spectrum);
}

double R_PE_Migdal(unsigned int nPE, double mu_PE, double sigma_PE, const DM_Particle& DM, DM_Distribution& DM_distr, Element& element, Atom& atom, std::vector<double> electron_spectrum)
{
	if(electron_spectrum.empty())
		for(unsigned ne = 1; ne < 16; ne++)
			electron_spectrum.push_back(R_ne_Migdal(ne, DM, DM_distr, element, atom));
	return R_PE_Migdal_aux(nPE, mu_PE, sigma_PE, electron_spectrum);
}

//2. Electron recoil direct detection experiment with isolated target atoms
DM_Detector_Migdal::DM_Detector_Migdal()
: DM_Detector("Migdal scattering experiment", kg * year, "Nuclei"), target_element(Get_Element("Xe")), target_atom(Import_Ionization_Form_Factors("Xenon")), ne_threshold(0), ne_max(0), using_electron_threshold(false), using_electron_bins(false), PE_threshold(0), PE_max(0), S2_mu(0.0), S2_sigma(0.0), using_S2_threshold(false), using_S2_bins(false)
{
}

DM_Detector_Migdal::DM_Detector_Migdal(std::string label, double expo, std::string atom)
: DM_Detector(label, expo, "Nuclei"), target_element(Get_Element(atom)), target_atom(Import_Ionization_Form_Factors(atom)), ne_threshold(0), ne_max(0), using_electron_threshold(false), using_electron_bins(false), PE_threshold(0), PE_max(0), S2_mu(0.0), S2_sigma(0.0), using_S2_threshold(false), using_S2_bins(false)
{
}

//DM functions
double DM_Detector_Migdal::Maximum_Energy_Deposit(const DM_Particle& DM, const DM_Distribution& DM_distr) const
{
	double vMax = DM_distr.Maximum_DM_Speed();
	double Emax = 0.0;
	for(unsigned int j = 0; j < target_element.Number_of_Isotopes(); j++)
		for(auto& isotope : target_element.isotopes)
		{
			double ERmax = Maximum_Nuclear_Recoil_Energy(vMax, DM.mass, isotope.mass);
			if(ERmax > Emax && DM.Sigma_Nucleus(isotope, vMax) > 0.0)
				Emax = ERmax;
		}
	return Emax;
}

double DM_Detector_Migdal::Minimum_DM_Mass(DM_Particle& DM, const DM_Distribution& DM_distr) const
{
	double vMax	 = DM_distr.Maximum_DM_Speed();
	double E_min = target_atom.Lowest_Binding_Energy();
	if(using_electron_bins)
		E_min += (ne_threshold - 1.0) * target_atom.W;
	double mMin = 2.0 * E_min / vMax / vMax;
	return mMin;
}

double DM_Detector_Migdal::Minimum_DM_Speed(const DM_Particle& DM) const
{
	return sqrt(2.0 * target_atom.Lowest_Binding_Energy() / DM.mass);
}

double DM_Detector_Migdal::dRdE(double E, const DM_Particle& DM, DM_Distribution& DM_distr)
{
	double dRdE_tot = 0.0;
	for(auto& electron : target_atom.electrons)
		dRdE_tot += dRdEe_Migdal(E, DM, DM_distr, target_element, electron);
	return flat_efficiency * dRdE_tot;
}

double DM_Detector_Migdal::DM_Signals_Total(const DM_Particle& DM, DM_Distribution& DM_distr)
{
	double N = 0;

	if(statistical_analysis == "Binned Poisson")
	{
		std::vector<double> binned_events = DM_Signals_Binned(DM, DM_distr);
		N								  = std::accumulate(binned_events.begin(), binned_events.end(), 0.0);
	}
	else if(using_electron_threshold)
	{
		for(unsigned int ne = ne_threshold; ne <= ne_max; ne++)
		{
			N += flat_efficiency * exposure * R_ne_Migdal(ne, DM, DM_distr, target_element, target_atom);
		}
	}
	else if(using_energy_threshold)
	{

		for(unsigned int i = 0; i < target_atom.electrons.size(); i++)
		{
			double kMax = target_atom[i].k_max;
			double Emax = kMax * kMax / 2.0 / mElectron;
			if(Emax > energy_threshold)
			{
				std::function<double(double)> dNdE = [this, i, &DM, &DM_distr](double E) {
					return flat_efficiency * exposure * dRdEe_Migdal(E, DM, DM_distr, target_element, target_atom[i]);
				};
				double eps = libphysica::Find_Epsilon(dNdE, energy_threshold, Emax, 1e-6);
				N += libphysica::Integrate(dNdE, energy_threshold, Emax, eps);
			}
		}
	}
	else if(using_S2_threshold)
	{
		// Precompute the electron spectrum to speep up the computation of the S2 spectrum
		std::vector<double> electron_spectrum;
		for(unsigned int ne = 1; ne < 16; ne++)
			electron_spectrum.push_back(R_ne_Migdal(ne, DM, DM_distr, target_element, target_atom));
		for(unsigned int PE = PE_threshold; PE <= PE_max; PE++)
		{
			double PE_eff = 1.0;
			if(Trigger_Efficiency_PE.empty() == false)
				PE_eff *= Trigger_Efficiency_PE[PE - 1];
			if(Acceptance_Efficiency_PE.empty() == false)
				PE_eff *= Acceptance_Efficiency_PE[PE - 1];
			N += flat_efficiency * PE_eff * exposure * R_PE_Migdal(PE, S2_mu, S2_sigma, DM, DM_distr, target_element, target_atom, electron_spectrum);
		}
	}

	return N;
}

std::vector<double> DM_Detector_Migdal::DM_Signals_Binned(const DM_Particle& DM, DM_Distribution& DM_distr)
{
	if(statistical_analysis != "Binned Poisson")
	{
		std::cerr << "Error in obscura::DM_Detector_Ionization::DM_Signals_Binned(): Statistical analysis is " << statistical_analysis << ", not 'Binned Poisson'" << std::endl;
		std::exit(EXIT_FAILURE);
	}
	else if(using_energy_bins)
	{
		return DM_Signals_Energy_Bins(DM, DM_distr);
	}
	else if(using_electron_bins)
	{
		return DM_Signals_Electron_Bins(DM, DM_distr);
	}
	else if(using_S2_bins)
	{
		return DM_Signals_PE_Bins(DM, DM_distr);
	}
	else
	{
		std::cerr << "Error in obscura::DM_Detector_Ionization::DM_Signals_Binned(): Statistical analysis is 'Binned Poisson' but no bins have been defined. This should not happen ever." << std::endl;
		std::exit(EXIT_FAILURE);
	}
}

//Electron spectrum
std::vector<double> DM_Detector_Migdal::DM_Signals_Electron_Bins(const DM_Particle& DM, DM_Distribution& DM_distr)
{
	if(!using_electron_bins)
	{
		std::cerr << "Error in obscura::DM_Detector_Migdal::DM_Signals_Electron_Bins(const DM_Particle&,DM_Distribution&): Not using electron bins." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	else
	{
		std::vector<double> signals;
		for(unsigned int bin = 0; bin < number_of_bins; bin++)
		{
			unsigned int ne = ne_threshold + bin;
			double R_bin	= R_ne_Migdal(ne, DM, DM_distr, target_element, target_atom);
			signals.push_back(flat_efficiency * bin_efficiencies[bin] * exposure * R_bin);
		}
		return signals;
	}
}

void DM_Detector_Migdal::Use_Electron_Threshold(unsigned int ne_thr, unsigned int nemax)
{
	Initialize_Poisson();
	using_electron_threshold = true;
	using_S2_threshold		 = false;
	using_energy_threshold	 = false;

	ne_threshold = ne_thr;
	ne_max		 = (nemax > 0 && nemax > ne_thr) ? nemax : 15;
	if(ne_max < ne_threshold)
	{
		std::cerr << "Error in obscura::DM_Detector::Use_Electron_Threshold(): ne threshold (" << ne_threshold << ") is higher than maximum (" << ne_max << ")." << std::endl;
		std::exit(EXIT_FAILURE);
	}
}

void DM_Detector_Migdal::Use_Electron_Bins(unsigned int ne_thr, unsigned int N_bins)
{
	Initialize_Binned_Poisson(N_bins);
	using_electron_bins = true;

	ne_threshold = ne_thr;
	ne_max		 = ne_threshold + N_bins - 1;
	if(ne_max < ne_threshold)
	{
		std::cerr << "Error in obscura::DM_Detector::Use_Electron_Bins(): ne threshold (" << ne_threshold << ") is higher than maximum (" << ne_max << ")." << std::endl;
		std::exit(EXIT_FAILURE);
	}
}

//PE (or S2) spectrum
std::vector<double> DM_Detector_Migdal::DM_Signals_PE_Bins(const DM_Particle& DM, DM_Distribution& DM_distr)
{
	if(!using_S2_bins)
	{
		std::cerr << "Error in obscura::DM_Detector_Migdal::DM_Signals_PE_Bins(const DM_Particle&,DM_Distribution&): Not using PE bins." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	else
	{
		// Precompute the electron spectrum to speep up the computation of the S2 spectrum
		std::vector<double> electron_spectrum;
		for(unsigned int ne = 1; ne < 16; ne++)
			electron_spectrum.push_back(R_ne_Migdal(ne, DM, DM_distr, target_element, target_atom));

		std::vector<double> signals;
		for(unsigned int bin = 0; bin < number_of_bins; bin++)
		{
			double R_bin = 0.0;
			for(unsigned int nPE = S2_bin_ranges[bin]; nPE < S2_bin_ranges[bin + 1]; nPE++)
			{
				double R_new = R_PE_Migdal(nPE, S2_mu, S2_sigma, DM, DM_distr, target_element, target_atom, electron_spectrum);
				if(Trigger_Efficiency_PE.empty() == false)
					R_new *= Trigger_Efficiency_PE[nPE - 1];
				if(Acceptance_Efficiency_PE.empty() == false)
					R_new *= Acceptance_Efficiency_PE[nPE - 1];
				R_bin += R_new;
			}
			signals.push_back(flat_efficiency * bin_efficiencies[bin] * exposure * R_bin);
		}
		return signals;
	}
}

void DM_Detector_Migdal::Use_PE_Threshold(double S2mu, double S2sigma, unsigned int nPE_thr, unsigned int nPE_max)
{
	Initialize_Poisson();
	using_S2_threshold		 = true;
	using_electron_threshold = false;
	using_energy_threshold	 = false;

	S2_mu	 = S2mu;
	S2_sigma = S2sigma;

	PE_threshold = nPE_thr;
	PE_max		 = nPE_max;
	if(PE_max < PE_threshold)
	{
		std::cerr << "Error in obscura::DM_Detector::Use_PE_Threshold(): PE threshold (" << PE_threshold << ") is higher than maximum (" << PE_max << ")." << std::endl;
		std::exit(EXIT_FAILURE);
	}
}

void DM_Detector_Migdal::Import_Trigger_Efficiency_PE(std::string filename)
{
	if(!using_S2_bins && !using_S2_threshold)
	{
		std::cerr << "Error in obscura::DM_Detector_Migdal::Import_Trigger_Efficiency_PE(): No PE spectrum has been initialized." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	else
	{
		Trigger_Efficiency_PE = libphysica::Import_List(filename);
	}
}

void DM_Detector_Migdal::Import_Acceptance_Efficiency_PE(std::string filename)
{
	if(!using_S2_bins && !using_S2_threshold)
	{
		std::cerr << "Error in obscura::DM_Detector_Migdal::Import_Acceptance_Efficiency_PE(): No PE spectrum has been initialized." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	else
	{
		Acceptance_Efficiency_PE = libphysica::Import_List(filename);
	}
}

//Binned Poisson:  PE bins (S2)
void DM_Detector_Migdal::Use_PE_Bins(double S2mu, double S2sigma, const std::vector<unsigned int>& bin_ranges)
{
	Initialize_Binned_Poisson(bin_ranges.size() - 1);
	using_S2_bins = true;

	S2_mu		  = S2mu;
	S2_sigma	  = S2sigma;
	S2_bin_ranges = bin_ranges;
}

void DM_Detector_Migdal::Print_Summary(int MPI_rank) const
{
	Print_Summary_Base();
	std::cout << std::endl
			  << "\tMigdal scattering experiment." << std::endl
			  << "\tTarget:\t\t\t" << target_atom.name << std::endl
			  << "\tElectron bins:\t\t" << (using_electron_bins ? "[x]" : "[ ]") << std::endl
			  << "\tPE (S2) bins:\t\t" << (using_S2_bins ? "[x]" : "[ ]") << std::endl;
	if(using_S2_bins || using_S2_threshold)
	{
		std::cout << "\tmu_PE:\t\t" << S2_mu << std::endl
				  << "\tsigma_PE:\t" << S2_sigma << std::endl
				  << "\tImported trigger efficiencies:\t" << (Trigger_Efficiency_PE.empty() ? "[ ]" : "[x]") << std::endl
				  << "\tImported acc. efficiencies:\t" << (Acceptance_Efficiency_PE.empty() ? "[ ]" : "[x]") << std::endl;
		if(using_S2_bins)
		{
			std::cout << "\n\t\tBin\tBin range [S2]" << std::endl;
			for(unsigned int bin = 0; bin < number_of_bins; bin++)
			{
				std::cout << "\t\t" << bin + 1 << "\t[" << S2_bin_ranges[bin] << "," << S2_bin_ranges[bin + 1] << ")" << std::endl;
			}
		}
	}
	else if(using_electron_bins || using_electron_threshold)
	{
		std::cout << "\t\tNe threshold:\t" << ne_threshold << std::endl;
		std::cout << "\t\tNe max:\t\t" << ne_max << std::endl;
	}
	std::cout << "----------------------------------------" << std::endl
			  << std::endl;
}

}	// namespace obscura
