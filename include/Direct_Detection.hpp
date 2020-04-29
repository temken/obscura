#ifndef __Direct_Detection_hpp_
#define __Direct_Detection_hpp_

#include <string>

//Headers from libphys library
#include "Natural_Units.hpp"
#include "Numerics.hpp"

#include "DM_Particle.hpp"
#include "DM_Distribution.hpp"

// DM Detector base class, which includes the statistical methods.
class DM_Detector
{
	protected:
		std::string targets;
		double exposure, flat_efficiency;
		double energy_threshold, energy_max;

		virtual double Maximum_Energy_Deposit(const DM_Particle& DM, const DM_Distribution& DM_distr) const { return 0.0; };
		virtual double Minimum_DM_Mass(DM_Particle& DM, const DM_Distribution& DM_distr) const {return 0.0;};
		
		//Statistics
		std::string statistical_analysis;

		//a) Poisson statistics
		unsigned long int observed_events;
		double expected_background;

		//b) Binned Poisson statistics
		void Use_Binned_Poisson(unsigned bins);
		unsigned int number_of_bins;
		std::vector<double> bin_energies;
		std::vector<double> bin_efficiencies;
		std::vector<unsigned long int> bin_observed_events;
		std::vector<double> bin_expected_background;

		//c) Maximum gap a'la Yellin
		std::vector<double> maximum_gap_energy_data;
		double P_Value_Maximum_Gap(const DM_Particle& DM, DM_Distribution& DM_distr);

		void Print_Summary_Base(int MPI_rank = 0) const;
		
	public:
		std::string name;
		DM_Detector() :  targets("base targets"), exposure(0.0), flat_efficiency(1.0), energy_threshold(0), energy_max(0), statistical_analysis("Poisson"), observed_events(0), expected_background(0.0), number_of_bins(0), name("base name") {};
		DM_Detector(std::string label, double expo,std::string target_type) : targets(target_type), exposure(expo) , flat_efficiency(1.0), energy_threshold(0), energy_max(0), statistical_analysis("Poisson"), observed_events(0), expected_background(0.0), number_of_bins(0), name(label) {};

		void Set_Flat_Efficiency(double eff);

		//DM functions
		virtual double Minimum_DM_Speed(const DM_Particle& DM) const {return 0.0;};
		virtual double dRdE(double E, const DM_Particle& DM, DM_Distribution& DM_distr) { return 0.0;};
		virtual double DM_Signals_Total(const DM_Particle& DM, DM_Distribution& DM_distr);
		virtual std::vector<double> DM_Signals_Binned(const DM_Particle& DM, DM_Distribution& DM_distr);

		//Statistics
		double Log_Likelihood(const DM_Particle& DM, DM_Distribution& DM_distr);
		double Likelihood(const DM_Particle& DM, DM_Distribution& DM_distr);
		double P_Value(const DM_Particle& DM, DM_Distribution& DM_distr);
		
		//a) Poisson
		void Use_Poisson_Statistics();
		void Set_Observed_Events(unsigned long int N);
		void Set_Expected_Background(double B);
		
		//b) Binned Poisson
		void Use_Energy_Bins(double Emin, double Emax, int bins);
		void Set_Observed_Events(std::vector<unsigned long int> Ni);
		void Set_Bin_Efficiencies(const std::vector<double>& eff);
		void Set_Expected_Background(const std::vector<double>& Bi);
		
		//c) Maximum gap
		void Use_Maximum_Gap(std::string filename_energy_data,double dim = keV);

		//Limits/Constraints
		double Upper_Limit(DM_Particle& DM, DM_Distribution& DM_distr, double certainty = 0.95);
		std::vector<std::vector<double>> Upper_Limit_Curve(DM_Particle& DM, DM_Distribution& DM_distr, std::vector<double> masses, double certainty = 0.95);

		virtual void Print_Summary(int MPI_rank = 0) const {Print_Summary_Base(MPI_rank);};
};

#endif