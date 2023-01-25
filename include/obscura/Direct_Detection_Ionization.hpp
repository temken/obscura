#ifndef __Direct_Detection_Ionization_hpp_
#define __Direct_Detection_Ionization_hpp_

#include "obscura/DM_Distribution.hpp"
#include "obscura/DM_Particle.hpp"
#include "obscura/Direct_Detection.hpp"
#include "obscura/Target_Atom.hpp"

namespace obscura
{

class DM_Detector_Ionization : public DM_Detector
{
  private:
	std::vector<Atom> atomic_targets;
	std::vector<double> relative_mass_fractions;

	double Energy_Gap() const;
	double Lowest_W() const;

	// Electron spectrum
	unsigned int ne_threshold, ne_max;
	// (a) Poisson: Electron threshold
	bool using_electron_threshold;
	// (b) Binned Poisson: Electron bins
	bool using_electron_bins;
	std::vector<double> DM_Signals_Electron_Bins(const DM_Particle& DM, DM_Distribution& DM_distr);

	// PE (or S2) spectrum
	unsigned int PE_threshold, PE_max;
	// (a) Poisson: PE threshold (S2)
	bool using_S2_threshold;
	// (b) Binned Poisson: PE bins (S2)
	bool using_S2_bins;
	std::vector<unsigned int> S2_bin_ranges;
	std::vector<double> DM_Signals_PE_Bins(const DM_Particle& DM, DM_Distribution& DM_distr);

	// Computation of the S2 spectrum
	std::string S2_spectrum_method;

	// (1) "Poisson+Gauss": Use Poission and normal distribution to calculate the probability of observing a given number of PE (S2) following Essig et al.
	// And use tabulated trigger and acceptance efficiencies
	double S2_mu, S2_sigma;
	std::vector<double> Trigger_Efficiency_PE;
	std::vector<double> Acceptance_Efficiency_PE;
	double R_S2_Bin(unsigned int S2_1, unsigned int S2_2, const DM_Particle& DM, DM_Distribution& DM_distr, std::vector<double> electron_spectrum = {});

	// (2) "Response matrix": Use a response matrix to calculate the probability of observing a given number of PE (S2) following the method of the response matrix
	libphysica::Matrix response_matrix;
	std::vector<std::vector<double>> energy_ranges, s2_bin_info;
	std::vector<double> Compute_S2_Spectrum(const DM_Particle& DM, DM_Distribution& DM_distr);

  public:
	DM_Detector_Ionization(std::string label, double expo, std::string target_particles, std::string atom);
	DM_Detector_Ionization(std::string label, double expo, std::string target_particles, std::vector<std::string> atoms, std::vector<double> mass_fractions = {});

	// DM functions from the base class
	virtual double Maximum_Energy_Deposit(DM_Particle& DM, const DM_Distribution& DM_distr) const override;
	virtual double Minimum_DM_Speed(DM_Particle& DM) const override;
	virtual double Minimum_DM_Mass(DM_Particle& DM, const DM_Distribution& DM_distr) const override;

	virtual double dRdE(double E, const DM_Particle& DM, DM_Distribution& DM_distr) override;
	virtual double DM_Signals_Total(const DM_Particle& DM, DM_Distribution& DM_distr) override;
	virtual std::vector<double> DM_Signals_Binned(const DM_Particle& DM, DM_Distribution& DM_distr) override;

	// Energy spectrum
	virtual double dRdE_Ionization(double E, const DM_Particle& DM, DM_Distribution& DM_distr, const Nucleus& nucleus, Atomic_Electron& shell);
	double dRdE_Ionization(double E, const DM_Particle& DM, DM_Distribution& DM_distr, Atom& atom);

	// Electron spectrum
	double R_ne(unsigned int ne, const DM_Particle& DM, DM_Distribution& DM_distr, double W, const Nucleus& nucleus, Atomic_Electron& shell);
	double R_ne(unsigned int ne, const DM_Particle& DM, DM_Distribution& DM_distr, Atom& atom);
	double R_ne(unsigned int ne, const DM_Particle& DM, DM_Distribution& DM_distr);
	// (a) Poisson: Electron threshold
	void Use_Electron_Threshold(unsigned int ne_thr, unsigned int nemax = 0);
	// (b) Binned Poisson: Electron bins
	void Use_Electron_Bins(unsigned int ne_thr, unsigned int N_bins);

	// PE (or S2) spectrum
	void Initialize_S2_Spectrum(std::string method, double S2mu = 0.0, double S2sigma = 0.0);

	double R_S2(unsigned int S2, const DM_Particle& DM, DM_Distribution& DM_distr, double W, const Nucleus& nucleus, Atomic_Electron& shell, std::vector<double> electron_spectrum = {});
	double R_S2(unsigned int S2, const DM_Particle& DM, DM_Distribution& DM_distr, Atom& atom, std::vector<double> electron_spectrum = {});
	double R_S2(unsigned int S2, const DM_Particle& DM, DM_Distribution& DM_distr, std::vector<double> electron_spectrum = {});

	// (a) Poisson: PE threshold (S2)
	void Use_PE_Threshold(unsigned int nPE_thr, unsigned int nPE_max);
	void Import_Trigger_Efficiency_PE(std::string filename);
	void Import_Acceptance_Efficiency_PE(std::string filename);
	// (b) Binned Poisson: PE bins (S2)
	void Use_PE_Bins(const std::vector<unsigned int>& bin_ranges);

	virtual void Print_Summary(int MPI_rank = 0) const override;
};

}	// namespace obscura

#endif