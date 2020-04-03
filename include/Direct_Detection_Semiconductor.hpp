#ifndef __Direct_Detection_Semiconductor_hpp_
#define __Direct_Detection_Semiconductor_hpp_

#include <string>
#include <vector>

#include "Direct_Detection.hpp"
#include "DM_Distribution.hpp"
#include "DM_Particle.hpp"
#include "Target_Electron.hpp"


//1. Event spectra and rates
	extern double dRdEe_Semiconductor(double Ee, const DM_Particle& DM,DM_Distribution& DM_distr, const Semiconductor& target_crystal);
	extern double dRdQ_Semiconductor(int Q, const DM_Particle& DM,DM_Distribution& DM_distr, const Semiconductor& target_crystal);
	extern double R_total_Semiconductor(int Qthreshold, const DM_Particle& DM,DM_Distribution& DM_distr, const Semiconductor& target_crystal);

//2. Electron recoil direct detection experiment with semiconductor target
	class DM_Detector_Semiconductor : public DM_Detector
	{
		private:
			Semiconductor semiconductor_target;
			unsigned int Q_threshold;

			virtual double Maximum_Energy_Deposit(const DM_Particle& DM, const DM_Distribution& DM_distr) const override;
			virtual double Minimum_DM_Mass(DM_Particle& DM, const DM_Distribution& DM_distr) const override;
		
		public:
			DM_Detector_Semiconductor(std::string label, double expo, std::string crys, unsigned int Q_thr);

			virtual double dRdE(double E, const DM_Particle& DM, DM_Distribution& DM_distr) override;
			virtual double Minimum_DM_Speed(const DM_Particle& DM) const override;

			virtual double DM_Signals_Total(const DM_Particle& DM, DM_Distribution& DM_distr) override;

			void Use_Q_Bins(unsigned int Q_thr, unsigned int N_bins = 0);
			virtual std::vector<double> DM_Signals_Binned(const DM_Particle& DM, DM_Distribution& DM_distr) override;

			virtual void Print_Summary(int MPI_rank = 0) const override;
	};

#endif