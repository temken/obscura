#ifndef __Direct_Detection_hpp_
#define __Direct_Detection_hpp_

#include <string>

//Headers from libphys library
#include "Numerics.hpp"

#include "DM_Particle.hpp"
#include "DM_Distribution.hpp"

//1. Detector base class
class Detector
{
	protected:
		
		std::string name, targets;
		double exposure, flat_efficiency;
		unsigned long int observed_signals;

		virtual double Maximum_Energy_Deposit(const DM_Particle& DM, const DM_Distribution& DM_distr) const { return 0.0; };
		virtual double Minimum_DM_Mass(DM_Particle& DM, const DM_Distribution& DM_distr) const {return 0.0;};
		// Interpolation Spectrum_Base(const DM_Particle& DM, const DM_Distribution& DM_distr, double Emin,double Emax,unsigned int points = 100);
		void Print_Summary_Base() const;
		
	public:
		Detector() : name("Generic"), targets("default"), exposure(0.0), flat_efficiency(1.0), observed_signals(0) {};
		Detector(double expo,std::string target_type) : name("Generic"), targets("default"), exposure(expo) , flat_efficiency(1.0), observed_signals(0) {};

		void Set_Name(std::string n);
		void Set_Flat_Efficiency(double eff);
		void Set_Observed_Signals(unsigned long int n);

		virtual double dRdE(double E, const DM_Particle& DM, const DM_Distribution& DM_distr, double vDM = 1e-3) { return 0.0;};
		// virtual Interpolation Spectrum(const DM_Particle& DM, const DM_Distribution& DM_distr, unsigned int points = 200) {return Interpolation();};
		
		virtual double N_Signals(const DM_Particle& DM, const DM_Distribution& DM_distr) { return 0.0; };
		virtual double Likelihood(const DM_Particle& DM, const DM_Distribution& DM_distr) { return 0.0;};

		std::vector<std::vector<double>> Limit_Curve(DM_Particle& DM, const DM_Distribution& DM_distr, double mMin,double mMax, int points = 50,double certainty = 0.95);
		virtual double Upper_Bound(DM_Particle& DM, const DM_Distribution& DM_distr, double certainty = 0.95) {return 0.0;};

		virtual double Minimum_DM_Speed(double mDM) const {return 0.0;};
		
		virtual void Print_Summary() const {Print_Summary_Base();};
};

#endif