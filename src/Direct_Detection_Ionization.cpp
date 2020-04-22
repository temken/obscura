#include "Direct_Detection_Ionization.hpp"

#include <cmath>

//Headers from libphys library
#include "Statistics.hpp"
#include "Utilities.hpp"

//1. Event spectra and rates
	double dRdEe_Ionization(double Ee, const DM_Particle& DM, DM_Distribution& DM_distr, const Atomic_Electron& shell)
	{
		double prefactor = 1.0/shell.nucleus_mass * DM_distr.DM_density / DM.mass / Ee / 2.0;
		double vDM = 1.0e-3; 	//cancels in the product with dSigma_dq^2 
								//-> THIS NEEDS TO BE UPDATED WHEN IMPLEMENTING VELOCITY DEPENDING CROSS SECTIONS
		
		// // Integration over q
		// // (a) Using numerical integration function.
		// double vMax = DM_distr.Maximum_DM_Speed();
		// double qMin = DM.mass * vMax - sqrt(DM.mass * DM.mass * vMax * vMax - 2.0 * DM.mass * shell.binding_energy);
		// double qMax = DM.mass * vMax + sqrt(DM.mass * DM.mass * vMax * vMax - 2.0 * DM.mass * shell.binding_energy);
		// if(qMax > shell.q_max) qMax = shell.q_max; 

		// std::function<double(double)> integrand = [Ee,&DM, vDM, &DM_distr, &shell] (double q)
		// {
		// 	double vMin = vMinimal_Electrons(q, shell.binding_energy + Ee, DM.mass);
		// 	return q * DM.dSigma_dq2_Electron(q,vDM) * vDM * vDM * DM_distr.Eta_Function(vMin) * shell.Ionization_Form_Factor(q,Ee);
		// };
		// double eps = Find_Epsilon(integrand, qMin, qMax, 1.0e-4);
		// double integral = Integrate(integrand, qMin, qMax, eps);
		
		// (b) Simply summing up the table entries
		double integral = 0.0;
		double k = sqrt(2.0 * mElectron * Ee);
		int ki = std::round(log10(k / shell.k_min) / shell.dlogk);
		if(ki >= shell.Nk || ki < 0) 
		{
 			std::cerr <<"Warning in dRdEe_Ionization(double,const DM_Particle&,DM_Distribution&,const Atomic_Electron&): Index ki = "<<ki<<" out of bounds. Function returns 0.0."<<std::endl;
 			return 0.0;
		}
		double vMax = DM_distr.Maximum_DM_Speed();
		double qMin = DM.mass * vMax - sqrt(DM.mass * DM.mass * vMax * vMax - 2.0 * DM.mass * shell.binding_energy);
		double qMax = DM.mass * vMax + sqrt(DM.mass * DM.mass * vMax * vMax - 2.0 * DM.mass * shell.binding_energy);
		int qi_min = std::floor(log10(qMin / shell.q_min) / shell.dlogq);
		int qi_max = std::floor(log10(qMax / shell.q_min) / shell.dlogq);
		if(qi_max > shell.Nq) qi_max = shell.Nq;
		// int qi_min = 0;
		// int qi_max = shell.Nq;
		for(int qi = qi_min; qi < qi_max; qi++)
		{
			double q = shell.q_Grid[qi];
			double vMin = vMinimal_Electrons(q, shell.binding_energy + Ee, DM.mass);
			integral += log(10.0) * shell.dlogq * q * q * DM.dSigma_dq2_Electron(q,vDM) * vDM * vDM * DM_distr.Eta_Function(vMin) *shell.Ionization_Form_Factor(q,Ee);//* shell.Ionization_Form_Factor[ki][qi];//
		}

		return prefactor * integral;
	}

	double dRdEe_Ionization(double Ee, const DM_Particle& DM, DM_Distribution& DM_distr, const Atom& atom)
	{
		double result = 0.0;
		for(unsigned int i = 0; i < atom.electrons.size(); i++)
		{
			result += dRdEe_Ionization(Ee, DM, DM_distr, atom.electrons[i]);
		}
		return result;
	}
	
	double PDF_ne(unsigned int ne, double Ee, const Atomic_Electron& shell)
	{
		double fR=0.0;
		double NxNi = 0.2;
		double fe = (1.0-fR)/(1.0+NxNi);
		double neMax = shell.number_of_secondary_electrons + std::floor( Ee / shell.W);
		return PMF_Binomial(neMax,fe,ne-1);
	}

	double R_ne_Ionization(unsigned int ne, const DM_Particle& DM, DM_Distribution& DM_distr, const Atomic_Electron& shell)
	{
		double sum=0.0;
		for(unsigned int ki=0; ki < shell.Nk; ki++)
		{
			double k = shell.k_Grid[ki];
			double Ee = k*k / 2.0 / mElectron;
			sum += log(10.0) * shell.dlogk * k*k/mElectron * PDF_ne(ne,Ee,shell) * dRdEe_Ionization(Ee, DM, DM_distr, shell);
		}
		return sum;
	}

	double R_ne_Ionization(unsigned int ne, const DM_Particle& DM, DM_Distribution& DM_distr, const Atom& atom)
	{
		double result = 0.0;
		for(unsigned int i = 0; i < atom.electrons.size(); i++)
		{
			result += R_ne_Ionization(ne, DM, DM_distr, atom.electrons[i]);
		}
		return result;
	}

	double R_PE_Ionization(unsigned int nPE, double mu_PE, double sigma_PE, const DM_Particle& DM, DM_Distribution& DM_distr, const Atomic_Electron& shell)
	{
		double sum=0.0;
		for(int ne = 1; ne < 16; ne++)
		{
				sum += PDF_Gauss(nPE,mu_PE*ne,sqrt(ne)*sigma_PE) * R_ne_Ionization(ne, DM, DM_distr, shell);
		}
		return sum;
	}

	double R_PE_Ionization(unsigned int nPE, double mu_PE, double sigma_PE, const DM_Particle& DM, DM_Distribution& DM_distr, const Atom& atom)
	{
		double result = 0.0;
		for(unsigned int i = 0; i < atom.electrons.size(); i++)
		{
			result += R_PE_Ionization(nPE, mu_PE, sigma_PE, DM, DM_distr, atom.electrons[i]);
		}
		return result;
	}

//2. Electron recoil direct detection experiment with isolated target atoms
	
	double DM_Detector_Ionization::Maximum_Energy_Deposit(const DM_Particle& DM, const DM_Distribution& DM_distr) const
	{
		double vMax = DM_distr.Maximum_DM_Speed();
		return DM.mass / 2.0 * vMax * vMax;
	}

	double DM_Detector_Ionization::Minimum_DM_Mass(DM_Particle& DM, const DM_Distribution& DM_distr) const
	{
		double vMax = DM_distr.Maximum_DM_Speed();
		double E_min = target_atom.Lowest_Binding_Energy();
		if(using_electron_bins) E_min += (ne_threshold - 1.0) * target_atom.W;
		double mMin = 2.0 * E_min / vMax / vMax;
		return mMin;
	}

	DM_Detector_Ionization::DM_Detector_Ionization(std::string label, double expo, Atom& atom)
	: DM_Detector(label, expo, "Electrons"), target_atom(atom)
	{

	}

	double DM_Detector_Ionization::dRdE(double E, const DM_Particle& DM, DM_Distribution& DM_distr) 
	{
		return dRdEe_Ionization(E, DM, DM_distr, target_atom);
	}

	double DM_Detector_Ionization::Minimum_DM_Speed(const DM_Particle& DM) const
	{
		return sqrt(2.0 * target_atom.Lowest_Binding_Energy() / DM.mass);
	}

	double DM_Detector_Ionization::DM_Signals_Total(const DM_Particle& DM, DM_Distribution& DM_distr) 
	{
		double N=0;
		if(using_electron_bins)
		{
			for(int ne = ne_threshold ; ne < 16 ; ne++)
			{
				N += exposure * R_ne_Ionization(ne,DM, DM_distr, target_atom);
			}
		}
		else if(using_S2_bins)
		{
			for(int nPE = S2_bin_ranges.front(); nPE < S2_bin_ranges.back(); nPE++)
			{
				double N_PE = exposure * R_PE_Ionization(nPE,S2_mu, S2_sigma, DM, DM_distr, target_atom);
				if(Trigger_Efficiency_PE.empty() == false) N_PE *= Trigger_Efficiency_PE[nPE - 1];
				if(Acceptance_Efficiency_PE.empty() == false) N_PE *= Acceptance_Efficiency_PE[nPE - 1];
				N += N_PE;
			}
		}
		else
		{
			double R=0.0;
			for(unsigned int i = 0; i < target_atom.electrons.size(); i++)
			{
				double k_min = sqrt(2.0 * mElectron * energy_threshold);
				int ki_min = std::floor(log10(k_min / target_atom[i].k_min) / target_atom[i].dlogk);
				for(unsigned int ki = ki_min; ki < target_atom[i].Nk; ki++)
				{
					double k = target_atom[i].k_Grid[ki];
					double Ee = k*k / 2.0 / mElectron;
					R += log(10.0) * target_atom[i].dlogk * k*k/mElectron * dRdEe_Ionization(Ee, DM, DM_distr, target_atom[i]);
				}
			}
			N = exposure * R;
		}
		return flat_efficiency * N;
	}

	void DM_Detector_Ionization::Use_Electron_Bins(unsigned int ne_thr, unsigned int N_bins)
	{
		using_electron_bins = true;
		using_S2_bins = false;

		ne_threshold = ne_thr;
		number_of_bins = N_bins;
	}

	void DM_Detector_Ionization::Use_PE_Bins(double S2mu, double S2sigma, const std::vector<int> &bin_ranges)
	{
		using_S2_bins = true;
		using_electron_bins = false;

		S2_mu = S2mu;
		S2_sigma = S2sigma;
		S2_bin_ranges = bin_ranges;

		number_of_bins = S2_bin_ranges.size() - 1;
	}

	void DM_Detector_Ionization::Import_Trigger_Efficiency_PE(std::string filename)
	{
		Trigger_Efficiency_PE = Import_List(filename);
	}

	void DM_Detector_Ionization::Import_Acceptance_Efficiency_PE(std::string filename)
	{
		Acceptance_Efficiency_PE = Import_List(filename);
	}

	std::vector<double> DM_Detector_Ionization::DM_Signals_Binned(const DM_Particle& DM, DM_Distribution& DM_distr) 
	{
		std::vector<double> N_binned;
		for(unsigned int bin = 0 ; bin < number_of_bins ; bin++)
		{
			double R_bin = 0.0;
			if(using_electron_bins)
			{
				unsigned int ne = ne_threshold + bin;
				R_bin = R_ne_Ionization(ne, DM, DM_distr,target_atom);
			}
			else if(using_S2_bins)
			{
				for(unsigned int nPE = S2_bin_ranges[bin]; nPE < S2_bin_ranges[bin+1]; nPE++)
				{
					double R_new = R_PE_Ionization(nPE,S2_mu, S2_sigma, DM, DM_distr, target_atom);
					if(Trigger_Efficiency_PE.empty() == false) R_new *= Trigger_Efficiency_PE[nPE - 1];
					if(Acceptance_Efficiency_PE.empty() == false) R_new *= Acceptance_Efficiency_PE[nPE - 1];
					R_bin += R_new;
				}
			}
			else
			{
				double E_min = bin_energies[bin];
				double E_max = bin_energies[bin+1];
				std::function<double(double)> spectrum = [this, &DM, &DM_distr] (double E)
				{
					return dRdEe_Ionization(E, DM, DM_distr,target_atom);
				};
				double epsilon = Find_Epsilon(spectrum, E_min, E_max, 1e-4);
				R_bin =  Integrate(spectrum, E_min, E_max, epsilon);
			}
			if(bin_efficiencies.empty() == false) R_bin *= bin_efficiencies[bin];
			N_binned.push_back(flat_efficiency * exposure * R_bin);
		}
		return N_binned;
	}

	void DM_Detector_Ionization::Print_Summary(int MPI_rank) const 
	{
		Print_Summary_Base();
		std::cout 	<<std::endl<<"Electron scattering experiment."<<std::endl
					<<"Target:\t\t\t"	<<target_atom.name <<std::endl
					<<"Electron bins:\t\t" <<(using_electron_bins? "[x]" : "[ ]") <<std::endl
					<<"PE spectrum:\t\t" <<(using_S2_bins? "[x]" : "[ ]") <<std::endl;
		if(using_S2_bins)
		{
			std::cout 	<<"\tmu_PE:\t" <<S2_mu<<std::endl
			 			<<"\tsigma_PE:\t" <<S2_sigma<<std::endl
			 			<<"\tImported trigger efficiencies:\t" <<(Trigger_Efficiency_PE.empty()? "[ ]" : "[x]") <<std::endl
			 			<<"\tImported acceptance efficiencies:\t" <<(Acceptance_Efficiency_PE.empty()? "[ ]" : "[x]") <<std::endl;

			std::cout <<"\n\tBin\tBin range [S2]"<<std::endl;
			for(unsigned int bin = 0; bin < number_of_bins; bin++)
			{
				std::cout<<"\t"<<bin+1<<"\t["<<S2_bin_ranges[bin]<<","<<S2_bin_ranges[bin+1]<<")"<<std::endl;
			}
		}
		else if(using_electron_bins)
		{
			std::cout<<"\tNe threshold:\t"<<ne_threshold<<std::endl;
		}
	 	std::cout<<"----------------------------------------"<<std::endl<<std::endl;
	}