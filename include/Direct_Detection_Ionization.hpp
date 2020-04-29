#ifndef __Direct_Detection_Ionization_hpp_
#define __Direct_Detection_Ionization_hpp_

// #include "Direct_Detection.hpp"
// #include "DM_Particle.hpp"
// #include "DM_Distribution.hpp"
// #include "Target_Electron.hpp"


// //1. Event spectra and rates
// 	extern double dRdEe_Ionization(double Ee, const DM_Particle& DM, DM_Distribution& DM_distr, const Atomic_Electron& shell);
// 	extern double dRdEe_Ionization(double Ee, const DM_Particle& DM, DM_Distribution& DM_distr, const Atom& atom);
	
// 	extern double R_ne_Ionization(unsigned int ne, const DM_Particle& DM, DM_Distribution& DM_distr, const Atomic_Electron& shell);
// 	extern double R_ne_Ionization(unsigned int ne, const DM_Particle& DM, DM_Distribution& DM_distr, const Atom& atom);
	
// 	extern double R_PE_Ionization(unsigned int nPE, double mu_PE, double sigma_PE, const DM_Particle& DM, DM_Distribution& DM_distr, const Atomic_Electron& shell);
// 	extern double R_PE_Ionization(unsigned int nPE, double mu_PE, double sigma_PE, const DM_Particle& DM, DM_Distribution& DM_distr, const Atom& atom);

// //2. Electron recoil direct detection experiment with isolated target atoms
// 	class DM_Detector_Ionization : public DM_Detector
// 	{
// 		private:
// 			Atom target_atom;			

// 			//Poisson with electron threshold
// 			bool using_electron_threshold;
// 			unsigned int ne_threshold;

// 			//Binned Poisson: Electron bins (n_e)
// 			bool using_electron_bins;
// 			std::vector<double> DM_Signals_Electron_Bins(const DM_Particle& DM, DM_Distribution& DM_distr);
			
// 			//Binned Poisson:  PE bins (S2)
// 			bool using_S2_bins;
// 			double S2_mu, S2_sigma;
// 			std::vector<int> S2_bin_ranges;
// 			std::vector<double> Trigger_Efficiency_PE;
// 			std::vector<double> Acceptance_Efficiency_PE;
// 			std::vector<double> DM_Signals_PE_Bins(const DM_Particle& DM, DM_Distribution& DM_distr);
			
// 			//DM functions
// 			virtual double Maximum_Energy_Deposit(const DM_Particle& DM, const DM_Distribution& DM_distr) const override;
// 			virtual double Minimum_DM_Mass(DM_Particle& DM, const DM_Distribution& DM_distr) const override;	


// 		public:
// 			DM_Detector_Ionization(std::string label, double expo, Atom& atom);

// 			//DM functions
// 			virtual double Minimum_DM_Speed(const DM_Particle& DM) const override;
// 			virtual double dRdE(double E, const DM_Particle& DM, DM_Distribution& DM_distr) override;
// 			virtual double DM_Signals_Total(const DM_Particle& DM, DM_Distribution& DM_distr) override;
// 			virtual std::vector<double> DM_Signals_Binned(const DM_Particle& DM, DM_Distribution& DM_distr) override;

// 			//Poisson with electron threshold
// 			virtual void Use_Poisson_Statistics() override;
// 			void Set_Electron_Threshold(unsigned int ne_thr);

// 			//Binned Poisson: Electron bins (n_e)
// 			void Use_Electron_Bins(unsigned int ne_thr, unsigned int N_bins);
			
// 			//Binned Poisson:  PE bins (S2)
// 			void Use_PE_Bins(double S2mu, double S2sigma, const std::vector<int> &bin_ranges);
// 			void Import_Trigger_Efficiency_PE(std::string filename);
// 			void Import_Acceptance_Efficiency_PE(std::string filename);

// 			virtual void Print_Summary(int MPI_rank = 0) const override;
// 	};

#endif