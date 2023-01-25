#include "obscura/Direct_Detection_Ionization.hpp"

#include <cmath>

#include "libphysica/Integration.hpp"
#include "libphysica/Natural_Units.hpp"
#include "libphysica/Special_Functions.hpp"
#include "libphysica/Statistics.hpp"
#include "libphysica/Utilities.hpp"

namespace obscura
{
using namespace libphysica::natural_units;

double DM_Detector_Ionization::Energy_Gap() const
{
	double gap = 1e10;
	for(auto& atom : atomic_targets)
		if(atom.Lowest_Binding_Energy() < gap)
			gap = atom.Lowest_Binding_Energy();
	return gap;
}

double DM_Detector_Ionization::Lowest_W() const
{
	double W = 1e10;
	for(auto& atom : atomic_targets)
		if(atom.W < W)
			W = atom.W;
	return W;
}

double DM_Detector_Ionization::Maximum_Energy_Deposit(DM_Particle& DM, const DM_Distribution& DM_distr) const
{
	double vMax = DM_distr.Maximum_DM_Speed();
	return DM.mass / 2.0 * vMax * vMax;
}
// Electron spectrum
std::vector<double> DM_Detector_Ionization::DM_Signals_Electron_Bins(const DM_Particle& DM, DM_Distribution& DM_distr)
{
	if(!using_electron_bins)
	{
		std::cerr << "Error in obscura::DM_Detector_Ionization::DM_Signals_Electron_Bins(const DM_Particle&,DM_Distribution&): Not using electron bins." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	else
	{
		std::vector<double> signals;
		for(unsigned int bin = 0; bin < number_of_bins; bin++)
		{
			unsigned int ne = ne_threshold + bin;
			double R_bin	= R_ne(ne, DM, DM_distr);
			signals.push_back(bin_efficiencies[bin] * exposure * R_bin);
		}
		return signals;
	}
}

// PE (or S2) spectrum
double DM_Detector_Ionization::R_S2_Bin(unsigned int S2_1, unsigned int S2_2, const DM_Particle& DM, DM_Distribution& DM_distr, std::vector<double> electron_spectrum)
{
	if(S2_spectrum_method == "Poisson+Gauss")
	{
		double R = 0.0;
		// Precompute the electron spectrum to speep up the computation of the S2 spectrum
		if(electron_spectrum.empty())
			for(unsigned int ne = 1; ne < 16; ne++)
				electron_spectrum.push_back(R_ne(ne, DM, DM_distr));
		for(unsigned int PE = S2_1; PE <= S2_2; PE++)
		{
			double PE_eff = 1.0;
			if(Trigger_Efficiency_PE.empty() == false)
				PE_eff *= Trigger_Efficiency_PE[PE - 1];
			if(Acceptance_Efficiency_PE.empty() == false)
				PE_eff *= Acceptance_Efficiency_PE[PE - 1];
			R += PE_eff * R_S2(PE, DM, DM_distr, electron_spectrum);
		}
		return R;
	}
	else if(S2_spectrum_method == "Response matrix")
	{
		std::vector<double> RS2_list = Compute_S2_Spectrum(DM, DM_distr);
		double R					 = 0;
		for(unsigned int i = 0; i < s2_bin_info.size(); i++)
		{
			unsigned int S2_min = s2_bin_info[i][0];
			unsigned int S2_max = s2_bin_info[i][1];
			if((S2_min > S2_1 && S2_min < S2_2) || (S2_max > S2_1 && S2_max < S2_2))
				R += RS2_list[i];
		}
		return R;
	}
	else
	{
		std::cerr << "Error in obscura::DM_Detector_Ionization::R_S2_Bin(unsigned int, unsigned int, const DM_Particle&, DM_Distribution&, std::vector<double>): Unknown S2 spectrum method: " << S2_spectrum_method << std::endl;
		std::exit(EXIT_FAILURE);
	}
}

std::vector<double> DM_Detector_Ionization::Compute_S2_Spectrum(const DM_Particle& DM, DM_Distribution& DM_distr)
{
	// Define the integrand, given by the energy spectrum.
	std::vector<double> dRdE_list(energy_ranges.size(), 0.0);
	std::function<double(double)> drde = [this, &DM, &DM_distr](double e) {
		return dRdE(e, DM, DM_distr);
	};

	// Integrate over the energy ranges of the response matrix
	for(unsigned int i = 0; i < energy_ranges.size(); i++)
	{
		double E_min = energy_ranges[i][1];
		double E_max = energy_ranges[i][2];
		double dR	 = libphysica::Integrate(drde, E_min, E_max);
		if(dR > 0.0)
			dRdE_list[i] = dR;
		else
			break;
	}

	// Multiply the response matrix with the energy spectrum
	libphysica::Vector dRdE_vector(dRdE_list);
	libphysica::Vector RS2_vector = response_matrix * dRdE_vector;
	std::vector<double> RS2_list(RS2_vector.Size());

	// Convert the vector to a list
	for(unsigned int i = 0; i < RS2_vector.Size(); i++)
		RS2_list[i] = RS2_vector[i];
	return RS2_list;
}

std::vector<double> DM_Detector_Ionization::DM_Signals_PE_Bins(const DM_Particle& DM, DM_Distribution& DM_distr)
{
	if(!using_S2_bins)
	{
		std::cerr << "Error in obscura::DM_Detector_Ionization::DM_Signals_PE_Bins(const DM_Particle&,DM_Distribution&): Not using PE bins." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	else if(S2_spectrum_method == "Poisson+Gauss")
	{
		// Precompute the electron spectrum to speep up the computation of the S2 spectrum
		std::vector<double> electron_spectrum;
		for(unsigned int ne = 1; ne < 16; ne++)
			electron_spectrum.push_back(R_ne(ne, DM, DM_distr));

		std::vector<double> signals;
		for(unsigned int bin = 0; bin < number_of_bins; bin++)
		{
			double R_bin = R_S2_Bin(S2_bin_ranges[bin], S2_bin_ranges[bin + 1] - 1, DM, DM_distr, electron_spectrum);
			signals.push_back(bin_efficiencies[bin] * exposure * R_bin);
		}
		return signals;
	}
	else if(S2_spectrum_method == "Response matrix")
	{
		std::vector<double> RS2_list = Compute_S2_Spectrum(DM, DM_distr);
		std::vector<double> signals(number_of_bins, 0.0);
		for(unsigned int i = 0; i < number_of_bins; i++)
		{
			double s2_min = S2_bin_ranges[i];
			double s2_max = S2_bin_ranges[i + 1];
			for(unsigned int j = 0; j < s2_bin_info.size(); j++)
				if(s2_bin_info[j][1] >= s2_min && s2_bin_info[j][2] < s2_max)
				{
					signals[i] += exposure * RS2_list[j];
					continue;
				}
		}
		return signals;
	}
	else
	{
		std::cerr << "Error in obscura::DM_Detector_Ionization::DM_Signals_PE_Bins(const DM_Particle&,DM_Distribution&): Unknown S2 spectrum method: " << S2_spectrum_method << std::endl;
		std::exit(EXIT_FAILURE);
	}
}

DM_Detector_Ionization::DM_Detector_Ionization(std::string label, double expo, std::string target_particle, std::string atom)
: DM_Detector(label, expo, target_particle), atomic_targets({atom}), relative_mass_fractions({1.0}), ne_threshold(3), ne_max(15), using_electron_threshold(false), using_electron_bins(false), PE_threshold(0), PE_max(0), using_S2_threshold(false), using_S2_bins(false), S2_spectrum_method("default"), S2_mu(0.0), S2_sigma(0.0)
{
}

DM_Detector_Ionization::DM_Detector_Ionization(std::string label, double expo, std::string target_particle, std::vector<std::string> atoms, std::vector<double> mass_fractions)
: DM_Detector(label, expo, target_particle), relative_mass_fractions(mass_fractions), ne_threshold(3), ne_max(15), using_electron_threshold(false), using_electron_bins(false), PE_threshold(0), PE_max(0), using_S2_threshold(false), using_S2_bins(false), S2_spectrum_method("default"), S2_mu(0.0), S2_sigma(0.0)

{
	for(auto& atom_name : atoms)
		atomic_targets.push_back(Atom(atom_name));
	if(relative_mass_fractions.size() == 0)
		relative_mass_fractions = std::vector<double>(atoms.size(), 1.0 / atoms.size());
	else if(relative_mass_fractions.size() != atomic_targets.size())
	{
		std::cerr << "Error in DM_Detector_Ionization::DM_Detector_Ionization(): Length of atomic_targets (" << atomic_targets.size() << ") and relative_mass_fractions (" << relative_mass_fractions.size() << ") does not match." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	else
	{
		double total = std::accumulate(relative_mass_fractions.begin(), relative_mass_fractions.end(), 0.0);
		if(total > 1.0)
			for(auto& entry : relative_mass_fractions)
				entry = entry / total;
	}
}

double DM_Detector_Ionization::Minimum_DM_Mass(DM_Particle& DM, const DM_Distribution& DM_distr) const
{
	double vMax = DM_distr.Maximum_DM_Speed();
	double E_min;
	if(using_energy_threshold)
		E_min = energy_threshold;
	else if(using_electron_bins || using_electron_threshold)
		E_min = Energy_Gap() + (ne_threshold - 1.0) * Lowest_W();
	else
		E_min = Energy_Gap();
	return 2.0 * E_min / vMax / vMax;
}

double DM_Detector_Ionization::Minimum_DM_Speed(DM_Particle& DM) const
{
	return sqrt(2.0 * Energy_Gap() / DM.mass);
}

double DM_Detector_Ionization::dRdE(double E, const DM_Particle& DM, DM_Distribution& DM_distr)
{
	double dRdE = 0.0;
	for(auto i = 0; i < atomic_targets.size(); i++)
		dRdE += relative_mass_fractions[i] * dRdE_Ionization(E, DM, DM_distr, atomic_targets[i]);
	return dRdE;
}

double DM_Detector_Ionization::DM_Signals_Total(const DM_Particle& DM, DM_Distribution& DM_distr)
{
	double N = 0;

	if(statistical_analysis == "Binned Poisson")
	{
		std::vector<double> binned_events = DM_Signals_Binned(DM, DM_distr);
		N								  = std::accumulate(binned_events.begin(), binned_events.end(), 0.0);
	}
	else if(using_electron_threshold)
		for(unsigned int ne = ne_threshold; ne <= ne_max; ne++)
			N += exposure * R_ne(ne, DM, DM_distr);
	else if(using_energy_threshold)
		for(auto i = 0; i < atomic_targets.size(); i++)
			for(auto& electron : atomic_targets[i].electrons)
			{
				double kMax = electron.k_max;
				double Emax = kMax * kMax / 2.0 / mElectron;
				if(Emax > energy_threshold)
				{
					std::function<double(double)> dNdE = [this, i, &electron, &DM, &DM_distr](double E) {
						return exposure * dRdE_Ionization(E, DM, DM_distr, atomic_targets[i].nucleus, electron);
					};
					N += relative_mass_fractions[i] * libphysica::Integrate(dNdE, energy_threshold, Emax);
				}
			}
	else if(using_S2_threshold)
		N = exposure * R_S2_Bin(PE_threshold, PE_max, DM, DM_distr);

	return N;
}

std::vector<double> DM_Detector_Ionization::DM_Signals_Binned(const DM_Particle& DM, DM_Distribution& DM_distr)
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

// Energy spectrum
double DM_Detector_Ionization::dRdE_Ionization(double E, const DM_Particle& DM, DM_Distribution& DM_distr, const Nucleus& nucleus, Atomic_Electron& shell)
{
	return 0.0;
}

double DM_Detector_Ionization::dRdE_Ionization(double E, const DM_Particle& DM, DM_Distribution& DM_distr, Atom& atom)
{
	double dRdE = 0.0;
	for(auto& electron : atom.electrons)
		dRdE += dRdE_Ionization(E, DM, DM_distr, atom.nucleus, electron);
	return dRdE;
}

// Electron spectrum

double PDF_ne(unsigned int ne, double Ee, double W, int n_secondary)
{
	double fR	 = 0.0;
	double NxNi	 = 0.2;
	double fe	 = (1.0 - fR) / (1.0 + NxNi);
	double neMax = n_secondary + std::floor(Ee / W);
	return libphysica::PMF_Binomial(neMax, fe, ne - 1);
}

double DM_Detector_Ionization::R_ne(unsigned int ne, const DM_Particle& DM, DM_Distribution& DM_distr, double W, const Nucleus& nucleus, Atomic_Electron& shell)
{
	double R = 0.0;
	for(auto& k : shell.k_Grid)
	{
		double Ee = k * k / 2.0 / mElectron;
		R += log(10.0) * shell.dlogk * k * k / mElectron * PDF_ne(ne, Ee, W, shell.number_of_secondary_electrons) * dRdE_Ionization(Ee, DM, DM_distr, nucleus, shell);
	}
	return R;
}

double DM_Detector_Ionization::R_ne(unsigned int ne, const DM_Particle& DM, DM_Distribution& DM_distr, Atom& atom)
{
	double R = 0.0;
	for(auto& electron : atom.electrons)
		R += R_ne(ne, DM, DM_distr, atom.W, atom.nucleus, electron);
	return R;
}

double DM_Detector_Ionization::R_ne(unsigned int ne, const DM_Particle& DM, DM_Distribution& DM_distr)
{
	double R = 0.0;
	for(auto i = 0; i < atomic_targets.size(); i++)
		R += relative_mass_fractions[i] * R_ne(ne, DM, DM_distr, atomic_targets[i]);
	return R;
}

void DM_Detector_Ionization::Use_Electron_Threshold(unsigned int ne_thr, unsigned int nemax)
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

void DM_Detector_Ionization::Use_Electron_Bins(unsigned int ne_thr, unsigned int N_bins)
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

// PE (or S2) spectrum
void DM_Detector_Ionization::Initialize_S2_Spectrum(std::string method, double S2mu, double S2sigma)
{
	S2_spectrum_method = method;
	if(S2_spectrum_method == "Poisson+Gauss")
	{
		S2_mu	 = S2mu;
		S2_sigma = S2sigma;
	}
	else if(S2_spectrum_method == "Response matrix")
	{
		// Import response matrix
		std::vector<std::vector<double>> data = libphysica::Import_Table(TOP_LEVEL_DIR "data/" + name + "/s2_response_er.dat");
		response_matrix						  = libphysica::Matrix(data);
		response_matrix						  = response_matrix.Transpose();
		// Import energy and S2 ranges of response matrix
		energy_ranges = libphysica::Import_Table(TOP_LEVEL_DIR "data/" + name + "/s2_response_er_energybins.dat", {keV, keV, keV});
		s2_bin_info	  = libphysica::Import_Table(TOP_LEVEL_DIR "data/" + name + "/s2_binning_info.dat", {}, 1);
	}
	else
	{
		std::cerr << "Error in obscura::DM_Detector_Ionization::Initialize_S2_Spectrum(): Unknown method '" << S2_spectrum_method << "'." << std::endl;
		std::exit(EXIT_FAILURE);
	}
}

double R_S2_aux(unsigned int nPE, double mu_PE, double sigma_PE, std::vector<double> R_ne_spectrum)
{
	double sum = 0.0;
	for(int ne = 1; ne < 16; ne++)
		sum += libphysica::PDF_Gauss(nPE, mu_PE * ne, sqrt(ne) * sigma_PE) * R_ne_spectrum[ne - 1];
	return sum;
}

double DM_Detector_Ionization::R_S2(unsigned int S2, const DM_Particle& DM, DM_Distribution& DM_distr, double W, const Nucleus& nucleus, Atomic_Electron& shell, std::vector<double> electron_spectrum)
{
	if(electron_spectrum.empty())
		for(unsigned ne = 1; ne < 16; ne++)
			electron_spectrum.push_back(R_ne(ne, DM, DM_distr, W, nucleus, shell));
	return R_S2_aux(S2, S2_mu, S2_sigma, electron_spectrum);
}

double DM_Detector_Ionization::R_S2(unsigned int S2, const DM_Particle& DM, DM_Distribution& DM_distr, Atom& atom, std::vector<double> electron_spectrum)
{
	if(electron_spectrum.empty())
		for(unsigned ne = 1; ne < 16; ne++)
			electron_spectrum.push_back(R_ne(ne, DM, DM_distr, atom));
	return R_S2_aux(S2, S2_mu, S2_sigma, electron_spectrum);
}

double DM_Detector_Ionization::R_S2(unsigned int S2, const DM_Particle& DM, DM_Distribution& DM_distr, std::vector<double> electron_spectrum)
{
	if(electron_spectrum.empty())
		for(unsigned ne = 1; ne < 16; ne++)
			electron_spectrum.push_back(R_ne(ne, DM, DM_distr));
	return R_S2_aux(S2, S2_mu, S2_sigma, electron_spectrum);
}

void DM_Detector_Ionization::Use_PE_Threshold(unsigned int nPE_thr, unsigned int nPE_max)
{
	if(S2_spectrum_method == "default")
	{
		std::cerr << "Error in obscura::DM_Detector_Ionization::Use_PE_Threshold(): No S2 spectrum defined." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	Initialize_Poisson();
	using_S2_threshold		 = true;
	using_electron_threshold = false;
	using_energy_threshold	 = false;

	PE_threshold = nPE_thr;
	PE_max		 = nPE_max;
	if(PE_max < PE_threshold)
	{
		std::cerr << "Error in obscura::DM_Detector::Use_PE_Threshold(): PE threshold (" << PE_threshold << ") is higher than maximum (" << PE_max << ")." << std::endl;
		std::exit(EXIT_FAILURE);
	}
}

void DM_Detector_Ionization::Import_Trigger_Efficiency_PE(std::string filename)
{
	if(!using_S2_bins && !using_S2_threshold)
	{
		std::cerr << "Error in obscura::DM_Detector_Ionization::Import_Trigger_Efficiency_PE(): No PE spectrum has been initialized." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	else if(S2_spectrum_method != "Poisson+Gauss")
	{
		std::cerr << "Error in obscura::DM_Detector_Ionization::Import_Trigger_Efficiency_PE(): Not supported for " << S2_spectrum_method << std::endl;
		std::exit(EXIT_FAILURE);
	}
	else
		Trigger_Efficiency_PE = libphysica::Import_List(filename);
}

void DM_Detector_Ionization::Import_Acceptance_Efficiency_PE(std::string filename)
{
	if(!using_S2_bins && !using_S2_threshold)
	{
		std::cerr << "Error in obscura::DM_Detector_Ionization::Import_Acceptance_Efficiency_PE(): No PE spectrum has been initialized." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	else if(S2_spectrum_method != "Poisson+Gauss")
	{
		std::cerr << "Error in obscura::DM_Detector_Ionization::Import_Acceptance_Efficiency_PE(): Not supported for " << S2_spectrum_method << std::endl;
		std::exit(EXIT_FAILURE);
	}
	else
		Acceptance_Efficiency_PE = libphysica::Import_List(filename);
}

// Binned Poisson:  PE bins (S2)
void DM_Detector_Ionization::Use_PE_Bins(const std::vector<unsigned int>& bin_ranges)
{
	if(S2_spectrum_method == "default")
	{
		std::cerr << "Error in obscura::DM_Detector_Ionization::Use_PE_Threshold(): No S2 spectrum defined." << std::endl;
		std::exit(EXIT_FAILURE);
	}

	Initialize_Binned_Poisson(bin_ranges.size() - 1);
	using_S2_bins = true;

	S2_bin_ranges = bin_ranges;
	PE_threshold  = bin_ranges[0];
	PE_max		  = bin_ranges[bin_ranges.size() - 1];
}

void DM_Detector_Ionization::Print_Summary(int MPI_rank) const
{
	Print_Summary_Base();
	std::cout << std::endl
			  << "\tElectron recoil experiment (ionization)." << std::endl
			  << "\tTarget(s):\t\t\t" << std::endl;
	for(auto i = 0; i < atomic_targets.size(); i++)
		std::cout << "\t\t\t" << atomic_targets[i].nucleus.name << "\t(" << libphysica::Round(100.0 * relative_mass_fractions[i]) << "%)" << std::endl;
	std::cout << "\tElectron bins:\t\t" << (using_electron_bins ? "[x]" : "[ ]") << std::endl
			  << "\tPE (S2) bins:\t\t" << (using_S2_bins ? "[x]" : "[ ]") << std::endl;
	if(using_S2_bins || using_S2_threshold)
	{
		std::cout << "S2 spectrum:\t\t\t" << S2_spectrum_method << std::endl;
		if(S2_spectrum_method == "Poisson+Gauss")
			std::cout << "\tmu_PE:\t\t" << S2_mu << std::endl
					  << "\tsigma_PE:\t" << S2_sigma << std::endl
					  << "\tImported trigger efficiencies:\t" << (Trigger_Efficiency_PE.empty() ? "[ ]" : "[x]") << std::endl
					  << "\tImported acc. efficiencies:\t" << (Acceptance_Efficiency_PE.empty() ? "[ ]" : "[x]") << std::endl;
		else if(S2_spectrum_method == "Response matrix")
		{
		}
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
