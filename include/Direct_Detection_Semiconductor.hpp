#ifndef __Direct_Detection_Semiconductor_hpp_
#define __Direct_Detection_Semiconductor_hpp_

#include <string>
#include <vector>

#include "Direct_Detection.hpp"
#include "DM_Distribution.hpp"
#include "DM_Particle.hpp"
#include "Target_Electron.hpp"

namespace obscura
{

//1. Event spectra and rates
	extern double dRdEe_Semiconductor(double Ee, const DM_Particle& DM,DM_Distribution& DM_distr, const Semiconductor& target_crystal);
	extern double R_Q_Semiconductor(int Q, const DM_Particle& DM,DM_Distribution& DM_distr, const Semiconductor& target_crystal);
	extern double R_total_Semiconductor(int Qthreshold, const DM_Particle& DM,DM_Distribution& DM_distr, const Semiconductor& target_crystal);

//2. Electron recoil direct detection experiment with semiconductor target
	class DM_Detector_Semiconductor : public DM_Detector
	{
		private:
			Semiconductor semiconductor_target;

		//DM functions
			virtual double Maximum_Energy_Deposit(const DM_Particle& DM, const DM_Distribution& DM_distr) const override;
			virtual double Minimum_DM_Mass(DM_Particle& DM, const DM_Distribution& DM_distr) const override;

		//Q spectrum
			unsigned int Q_threshold;
			
			// (a) Poisson: Energy threshold
			bool using_Q_threshold;
			
			// (b) Binned Poisson: Energy bins
			bool using_Q_bins;
			std::vector<double> DM_Signals_Q_Bins(const DM_Particle& DM, DM_Distribution& DM_distr);

		public:
			DM_Detector_Semiconductor(std::string label, double expo, std::string crys);

		//DM functions
			virtual double Minimum_DM_Speed(const DM_Particle& DM) const override;
			virtual double dRdE(double E, const DM_Particle& DM, DM_Distribution& DM_distr) override;
			virtual double DM_Signals_Total(const DM_Particle& DM, DM_Distribution& DM_distr) override;
			virtual std::vector<double> DM_Signals_Binned(const DM_Particle& DM, DM_Distribution& DM_distr) override;

		//Q spectrum
			// (a) Poisson
			void Use_Q_Threshold(unsigned int Q_thr);
			// (b) Binned Poisson
			void Use_Q_Bins(unsigned int Q_thr, unsigned int N_bins = 0);

			virtual void Print_Summary(int MPI_rank = 0) const override;
	};

}	// namespace obscura

#endif