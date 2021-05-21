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

	//DM functions
	virtual double Maximum_Energy_Deposit(const DM_Particle& DM, const DM_Distribution& DM_distr) const override;

	//Electron spectrum
	unsigned int ne_threshold, ne_max;
	// (a) Poisson: Electron threshold
	bool using_electron_threshold;
	// (b) Binned Poisson: Electron bins
	bool using_electron_bins;
	std::vector<double> DM_Signals_Electron_Bins(const DM_Particle& DM, DM_Distribution& DM_distr);

	//PE (or S2) spectrum
	unsigned int PE_threshold, PE_max;
	double S2_mu, S2_sigma;
	std::vector<double> Trigger_Efficiency_PE;
	std::vector<double> Acceptance_Efficiency_PE;
	// (a) Poisson: PE threshold (S2)
	bool using_S2_threshold;
	// (b) Binned Poisson: PE bins (S2)
	bool using_S2_bins;
	std::vector<unsigned int> S2_bin_ranges;
	double R_S2_Bin(unsigned int S2_1, unsigned int S2_2, const DM_Particle& DM, DM_Distribution& DM_distr, std::vector<double> electron_spectrum = {});
	std::vector<double> DM_Signals_PE_Bins(const DM_Particle& DM, DM_Distribution& DM_distr);

  public:
	DM_Detector_Ionization();
	DM_Detector_Ionization(std::string label, double expo, std::string target_particles, std::string atom);
	DM_Detector_Ionization(std::string label, double expo, std::string target_particles, std::vector<std::string> atoms, std::vector<double> mass_fractions = {});

	//DM functions from the base class
	virtual double Minimum_DM_Speed(const DM_Particle& DM) const override;
	virtual double Minimum_DM_Mass(DM_Particle& DM, const DM_Distribution& DM_distr) const override;

	virtual double dRdE(double E, const DM_Particle& DM, DM_Distribution& DM_distr) override;
	virtual double DM_Signals_Total(const DM_Particle& DM, DM_Distribution& DM_distr) override;
	virtual std::vector<double> DM_Signals_Binned(const DM_Particle& DM, DM_Distribution& DM_distr) override;

	//Energy spectrum
	virtual double dRdE_Ionization(double E, const DM_Particle& DM, DM_Distribution& DM_distr, const Nucleus& nucleus, Atomic_Electron& shell);
	double dRdE_Ionization(double E, const DM_Particle& DM, DM_Distribution& DM_distr, Atom& atom);

	//Electron spectrum
	double R_ne(unsigned int ne, const DM_Particle& DM, DM_Distribution& DM_distr, double W, const Nucleus& nucleus, Atomic_Electron& shell);
	double R_ne(unsigned int ne, const DM_Particle& DM, DM_Distribution& DM_distr, Atom& atom);
	double R_ne(unsigned int ne, const DM_Particle& DM, DM_Distribution& DM_distr);
	// (a) Poisson: Electron threshold
	void Use_Electron_Threshold(unsigned int ne_thr, unsigned int nemax = 0);
	// (b) Binned Poisson: Electron bins
	void Use_Electron_Bins(unsigned int ne_thr, unsigned int N_bins);

	//PE (or S2) spectrum
	double R_S2(unsigned int S2, const DM_Particle& DM, DM_Distribution& DM_distr, double W, const Nucleus& nucleus, Atomic_Electron& shell, std::vector<double> electron_spectrum = {});
	double R_S2(unsigned int S2, const DM_Particle& DM, DM_Distribution& DM_distr, Atom& atom, std::vector<double> electron_spectrum = {});
	double R_S2(unsigned int S2, const DM_Particle& DM, DM_Distribution& DM_distr, std::vector<double> electron_spectrum = {});

	// (a) Poisson: PE threshold (S2)
	void Use_PE_Threshold(double S2mu, double S2sigma, unsigned int nPE_thr, unsigned int nPE_max);
	void Import_Trigger_Efficiency_PE(std::string filename);
	void Import_Acceptance_Efficiency_PE(std::string filename);
	// (b) Binned Poisson: PE bins (S2)
	void Use_PE_Bins(double S2mu, double S2sigma, const std::vector<unsigned int>& bin_ranges);

	virtual void Print_Summary(int MPI_rank = 0) const override;
};

}	// namespace obscura

#endif