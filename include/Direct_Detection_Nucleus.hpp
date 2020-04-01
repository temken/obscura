#ifndef __Direct_Detection_Nucleus_hpp_
#define __Direct_Detection_Nucleus_hpp_

#include <vector>
#include <string>

//Headers from libphys library
#include "Natural_Units.hpp"
#include "Numerics.hpp"

#include "Direct_Detection.hpp"
#include "DM_Particle.hpp"
#include "DM_Distribution.hpp"

//1. Theoretical nuclear recoil spectrum [events per time, energy, and target mass]
	extern double dRdER_Nucleus(double ER, const DM_Particle& DM, DM_Distribution& DM_distr, const Isotope& target_isotope);
	extern double dRdER_Nucleus(double ER, const DM_Particle& DM, DM_Distribution& DM_distr, const Element& target_element);

//2. Nuclear recoil direct detection experiment
	class DM_Detector_Nucleus : public DM_Detector
	{
		private:
			//Target material
			std::vector<Element> target_elements;
			std::vector<double> relative_mass_fractions;

			//Experimental parameters
			double energy_resolution;
			bool using_efficiency_tables;
			std::vector<Interpolation> efficiencies;

			virtual double Maximum_Energy_Deposit(const DM_Particle& DM, const DM_Distribution& DM_distr) const override;
			virtual double Minimum_DM_Mass(DM_Particle& DM, const DM_Distribution& DM_distr) const override;
		
		public:
			DM_Detector_Nucleus();
			DM_Detector_Nucleus(std::string label, double expo,std::vector<Element> elements, double thr,double emax,std::vector<double> abund = {});

			void Set_Resolution(double res);
			void Import_Efficiency(std::string filename,double dim = keV);
			void Import_Efficiency(std::vector<std::string> filenames,double dim = keV);
			
			virtual double dRdE(double E, const DM_Particle& DM, DM_Distribution& DM_distr) override;
			virtual double Minimum_DM_Speed(const DM_Particle& DM) const override;
			
			virtual void Print_Summary(int MPI_rank = 0) const override;
	};
#endif