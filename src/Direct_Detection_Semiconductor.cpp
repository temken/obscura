#include "Direct_Detection_Semiconductor.hpp"

#include <cmath>

//1. Event spectra and rates

	double Minimum_Electron_Energy(int Q, const Semiconductor& target)
	{
		return target.epsilon * (Q - 1) + target.energy_gap;
	}

	unsigned int Electron_Hole_Pairs(double Ee, const Semiconductor& target)
	{
		return std::floor( (Ee-target.energy_gap) / target.epsilon +1 );
	}

	double dRdEe_Semiconductor(double Ee, const DM_Particle& DM,DM_Distribution& DM_distr, const Semiconductor& target_crystal)
	{
		//Integrate by summing over the tabulated form factors
		int Ei = std::round( Ee/target_crystal.dE -1 );
		if(Ei < 0 || Ei > 499)
		{
			std::cerr <<"Error in dRdEe_Semiconductor(): Ee lies beyond the tabulated cyrstal form factor."<<std::endl;
			std::exit(EXIT_FAILURE);
		}
		else
		{
			double prefactor = DM_distr.DM_density / DM.mass / target_crystal.M_cell * aEM * mElectron * mElectron;
			double vDM = 1e-3; //cancels in v^2 * dSigma/dq^2
			double sum=0.0;
			for(int qi=0; qi<900; qi++) 
			{
				double q = (qi+1) * target_crystal.dq;
				sum += target_crystal.dq / q / q * DM_distr.Eta_Function(vMinimal_Electrons(q,Ee,DM.mass)) * target_crystal.Crystal_Form_Factor[qi][Ei] * 4.0*vDM*vDM * DM.dSigma_dq2_Electron(q,vDM);
			}
			return prefactor * sum;
		}
		
	}

	double dRdQ_Semiconductor(int Q, const DM_Particle& DM,DM_Distribution& DM_distr, const Semiconductor& target_crystal)
	{
		//Energy threshold
		double Emin = Minimum_Electron_Energy(Q, target_crystal);
		double Emax = Minimum_Electron_Energy(Q+1, target_crystal);
		//Integrate over energies
		double sum =0.0;
		for(int Ei = (Emin/target_crystal.dE); Ei<500; Ei++)
		{
			double E = (Ei+1) * target_crystal.dE;
			if(E > Emax) break;
			// std::cout <<"\tAdd " <<Ei<<"\t"<<E/eV <<" between ["<<Emin/eV<<","<<Emax/eV<<"]"<<std::endl;
			sum+=target_crystal.dE * dRdEe_Semiconductor(E, DM, DM_distr, target_crystal);
		}
		return sum;
	}

	double R_total_Semiconductor(int Qthreshold, const DM_Particle& DM,DM_Distribution& DM_distr, const Semiconductor& target_crystal)
	{
		//Energy threshold
		double E_min = Minimum_Electron_Energy(Qthreshold, target_crystal);
		//Integrate over energies
		double sum =0.0;
		for(int Ei = (E_min /target_crystal.dE); Ei<500; Ei++)
		{
			double E = (Ei+1) * target_crystal.dE;
			sum+=target_crystal.dE * dRdEe_Semiconductor(E, DM, DM_distr, target_crystal);
		}
		return sum;
	}


//2. Electron recoil direct detection experiment with semiconductor target
	DM_Detector_Semiconductor::DM_Detector_Semiconductor(std::string label, double expo, std::string crys, unsigned int Q_thr)
	: DM_Detector(label, expo, "Electrons"), semiconductor_target(Semiconductor(crys)), Q_threshold(Q_thr)
	{
		energy_threshold= Minimum_Electron_Energy(Q_threshold, semiconductor_target);
		energy_max= Minimum_Electron_Energy(semiconductor_target.Q_max, semiconductor_target);
	}

	double DM_Detector_Semiconductor::Maximum_Energy_Deposit(const DM_Particle& DM, const DM_Distribution& DM_distr) const
	{
		return DM.mass/2.0 * pow(DM_distr.Maximum_DM_Speed(),2.0);
	}

	double DM_Detector_Semiconductor::Minimum_DM_Mass(DM_Particle& DM, const DM_Distribution& DM_distr) const
	{
		return 2.0 * energy_threshold * pow(DM_distr.Maximum_DM_Speed(), -2.0);
	}

	double DM_Detector_Semiconductor::dRdE(double E, const DM_Particle& DM, DM_Distribution& DM_distr)
	{
		return flat_efficiency * dRdEe_Semiconductor(E, DM, DM_distr, semiconductor_target);
	}

	double DM_Detector_Semiconductor::Minimum_DM_Speed(const DM_Particle& DM) const
	{
		return sqrt( 2.0 * energy_threshold / DM.mass);
	}

	double DM_Detector_Semiconductor::DM_Signals_Total(const DM_Particle& DM, DM_Distribution& DM_distr)
	{
		return exposure * flat_efficiency * R_total_Semiconductor(Q_threshold, DM, DM_distr, semiconductor_target);
	}

	void DM_Detector_Semiconductor::Use_Q_Bins(unsigned int Q_thr, unsigned int N_bins)
	{
		Q_threshold = Q_thr;
		unsigned int Q_max = (N_bins == 0)? semiconductor_target.Q_max : N_bins - Q_threshold + 1;
		Use_Binned_Poisson(Q_max - Q_threshold + 1);
		
		energy_threshold = Minimum_Electron_Energy(Q_threshold, semiconductor_target);
		energy_max = Minimum_Electron_Energy(Q_max, semiconductor_target);
	}

	std::vector<double> DM_Detector_Semiconductor::DM_Signals_Binned(const DM_Particle& DM, DM_Distribution& DM_distr)
	{
		std::vector<double> signals;
		for(int Q = Q_threshold; Q < Q_threshold + number_of_bins; Q++)
		{
			signals.push_back(exposure * flat_efficiency * bin_efficiencies[Q-1] * dRdQ_Semiconductor(Q, DM, DM_distr, semiconductor_target));
		}
		return signals;
	}

	void DM_Detector_Semiconductor::Print_Summary(int MPI_rank) const
	{
			Print_Summary_Base();
			std::cout 	<<std::endl <<"Electron scattering experiment."<<std::endl
						<<"Target:\t\t\t"	<<semiconductor_target.name <<" semiconductor"<<std::endl;
			std::cout 	<<"eh pair threshold:\t"<<Q_threshold<<std::endl
						<<"----------------------------------------"<<std::endl<<std::endl;
	}

