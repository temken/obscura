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
	DM_Detector_Ionization::DM_Detector_Ionization(std::string label, double expo, std::string atom)
	: DM_Detector(label, expo, "Electrons"), target_atom(Import_Ionization_Form_Factors(atom)),  ne_threshold(0), ne_max(0), using_electron_threshold(false), using_electron_bins(false), PE_threshold(0), PE_max(0), S2_mu (0.0), S2_sigma (0.0), using_S2_threshold(false), using_S2_bins(false)
	{
	}

	//DM functions
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

	double DM_Detector_Ionization::Minimum_DM_Speed(const DM_Particle& DM) const
	{
		return sqrt(2.0 * target_atom.Lowest_Binding_Energy() / DM.mass);
	}

	double DM_Detector_Ionization::dRdE(double E, const DM_Particle& DM, DM_Distribution& DM_distr) 
	{
		return flat_efficiency * dRdEe_Ionization(E, DM, DM_distr, target_atom);
	}


	double DM_Detector_Ionization::DM_Signals_Total(const DM_Particle& DM, DM_Distribution& DM_distr) 
	{
		double N = 0;
		if(statistical_analysis == "Binned Poisson")
		{
			std::vector<double> binned_events = DM_Signals_Binned(DM, DM_distr);
			N = std::accumulate(binned_events.begin(), binned_events.end(), 0.0);
		}
		else if(using_electron_threshold)
		{
			for(int ne = ne_threshold; ne <= ne_max; ne++)
			{
				N += flat_efficiency * exposure * R_ne_Ionization(ne, DM, DM_distr, target_atom);
			}
		}
		else if(using_energy_threshold)
		{
			for(unsigned int i = 0; i < target_atom.electrons.size(); i++)
			{
				double kMax = target_atom[i].k_max;
				double Emax = kMax*kMax/2.0/mElectron;
				if(Emax > energy_threshold)
				{
					std::function<double(double)> dNdE = [this, i, &DM, &DM_distr] (double E)
					{
						return flat_efficiency * exposure * dRdEe_Ionization(E, DM, DM_distr, target_atom[i]);
					};
					double eps =Find_Epsilon(dNdE,energy_threshold,Emax,1e-6);
					N +=  Integrate(dNdE,energy_threshold,Emax,eps);
				}
			}
		}
		else if(using_S2_threshold)
		{
			for(int PE = PE_threshold; PE <= PE_max; PE++)
			{
				double PE_eff = 1.0;
				if(Trigger_Efficiency_PE.empty() == false) PE_eff *= Trigger_Efficiency_PE[PE - 1];
				if(Acceptance_Efficiency_PE.empty() == false) PE_eff *= Acceptance_Efficiency_PE[PE - 1];
				N += flat_efficiency * PE_eff * exposure * R_PE_Ionization(PE, S2_mu, S2_sigma, DM, DM_distr, target_atom);
			}
		}
		return N;
	}

	std::vector<double> DM_Detector_Ionization::DM_Signals_Binned(const DM_Particle& DM, DM_Distribution& DM_distr) 
	{
		if(statistical_analysis != "Binned Poisson")
		{
			std::cerr <<"Error in DM_Detector_Ionization::DM_Signals_Binned(): Statistical analysis is "<<statistical_analysis <<", not 'Binned Poisson'" <<std::endl;
			std::exit(EXIT_FAILURE);
		}
		else if(using_energy_bins)
		{
			return DM_Signals_Energy_Bins(DM, DM_distr);
		}
		else if(using_electron_bins)
		{
			return DM_Signals_Electron_Bins(DM, DM_distr);
		}
		else if(using_S2_bins)
		{
			return DM_Signals_PE_Bins(DM, DM_distr);
		}
		else
		{
			std::cerr <<"Error in DM_Detector_Ionization::DM_Signals_Binned(): Statistical analysis is 'Binned Poisson' but no bins have been defined. This should not happen ever." <<std::endl;
			std::exit(EXIT_FAILURE);
		}
	}

	//Electron spectrum	
	std::vector<double> DM_Detector_Ionization::DM_Signals_Electron_Bins(const DM_Particle& DM, DM_Distribution& DM_distr)
	{
		if(!using_electron_bins)
		{
			std::cerr <<"Error in DM_Detector_Ionization::DM_Signals_Electron_Bins(const DM_Particle&,DM_Distribution&): Not using electron bins." <<std::endl;
			std::exit(EXIT_FAILURE);
		}
		else
		{
			std::vector<double> signals;
			for(unsigned int bin = 0; bin < number_of_bins; bin++)
			{
				unsigned int ne = ne_threshold + bin;
				double R_bin = R_ne_Ionization(ne, DM, DM_distr,target_atom);
				signals.push_back(flat_efficiency * bin_efficiencies[bin] * exposure * R_bin);
			}
			return signals;
		}
	}

	void DM_Detector_Ionization::Use_Electron_Threshold(unsigned int ne_thr,unsigned int nemax)
	{
		Initialize_Poisson();
		using_electron_threshold = true;
		using_S2_threshold = false;
		using_energy_threshold = false;

		ne_threshold = ne_thr;
		ne_max = (nemax > 0 && nemax > ne_thr)? nemax : 15;
		if(ne_max < ne_threshold)
		{
			std::cerr <<"Error in DM_Detector::Use_Electron_Threshold(): ne threshold (" <<ne_threshold <<") is higher than maximum (" <<ne_max<<")."<<std::endl;
			std::exit(EXIT_FAILURE);
		}
	}

	void DM_Detector_Ionization::Use_Electron_Bins(unsigned int ne_thr, unsigned int N_bins)
	{
		Initialize_Binned_Poisson(N_bins);
		using_electron_bins = true;
		
		ne_threshold = ne_thr;
		ne_max = ne_threshold + N_bins - 1;
		if(ne_max < ne_threshold)
		{
			std::cerr <<"Error in DM_Detector::Use_Electron_Bins(): ne threshold (" <<ne_threshold <<") is higher than maximum (" <<ne_max<<")."<<std::endl;
			std::exit(EXIT_FAILURE);
		}
	}

	//PE (or S2) spectrum
	std::vector<double> DM_Detector_Ionization::DM_Signals_PE_Bins(const DM_Particle& DM, DM_Distribution& DM_distr)
	{
		if(!using_S2_bins)
		{
			std::cerr <<"Error in DM_Detector_Ionization::DM_Signals_PE_Bins(const DM_Particle&,DM_Distribution&): Not using PE bins." <<std::endl;
			std::exit(EXIT_FAILURE);
		}
		else
		{
			std::vector<double> signals;
			for(unsigned int bin = 0; bin < number_of_bins; bin++)
			{
				double R_bin = 0.0;
				for(unsigned int nPE = S2_bin_ranges[bin]; nPE < S2_bin_ranges[bin+1]; nPE++)
				{
					double R_new = R_PE_Ionization(nPE, S2_mu, S2_sigma, DM, DM_distr, target_atom);
					if(Trigger_Efficiency_PE.empty() == false) R_new *= Trigger_Efficiency_PE[nPE - 1];
					if(Acceptance_Efficiency_PE.empty() == false) R_new *= Acceptance_Efficiency_PE[nPE - 1];
					R_bin += R_new;
				}				
				signals.push_back(flat_efficiency * bin_efficiencies[bin] * exposure * R_bin);
			}
			return signals;
		}
	}
	
	void DM_Detector_Ionization::Use_PE_Threshold(double S2mu, double S2sigma, unsigned int nPE_thr,unsigned int nPE_max)
	{
		Initialize_Poisson();
		using_S2_threshold = true;
		using_electron_threshold = false;
		using_energy_threshold = false;

		S2_mu = S2mu;
		S2_sigma = S2sigma;

		PE_threshold = nPE_thr;
		PE_max = nPE_max;
		if(PE_max < PE_threshold)
		{
			std::cerr <<"Error in DM_Detector::Use_PE_Threshold(): PE threshold (" <<PE_threshold <<") is higher than maximum (" <<PE_max<<")."<<std::endl;
			std::exit(EXIT_FAILURE);
		}
	}
	
	void DM_Detector_Ionization::Import_Trigger_Efficiency_PE(std::string filename)
	{	
		if(!using_S2_bins && !using_S2_threshold)
		{
			std::cerr <<"Error in DM_Detector_Ionization::Import_Trigger_Efficiency_PE(): No PE spectrum has been initialized." <<std::endl;
			std::exit(EXIT_FAILURE);
		}
		else
		{
			Trigger_Efficiency_PE = Import_List(filename);
		}
	}

	void DM_Detector_Ionization::Import_Acceptance_Efficiency_PE(std::string filename)
	{
		if(!using_S2_bins && !using_S2_threshold)
		{
			std::cerr <<"Error in DM_Detector_Ionization::Import_Acceptance_Efficiency_PE(): No PE spectrum has been initialized." <<std::endl;
			std::exit(EXIT_FAILURE);
		}
		else
		{
			Acceptance_Efficiency_PE = Import_List(filename);
		}
	}

	//Binned Poisson:  PE bins (S2)
	void DM_Detector_Ionization::Use_PE_Bins(double S2mu, double S2sigma, const std::vector<int> &bin_ranges)
	{
		Initialize_Binned_Poisson(bin_ranges.size() - 1);
		using_S2_bins = true;
		
		S2_mu = S2mu;
		S2_sigma = S2sigma;
		S2_bin_ranges = bin_ranges;
	}


	void DM_Detector_Ionization::Print_Summary(int MPI_rank) const 
	{
		Print_Summary_Base();
		std::cout 	<<std::endl<<"Electron scattering experiment."<<std::endl
					<<"Target:\t\t\t"	<<target_atom.name <<std::endl
					<<"Electron bins:\t\t" <<(using_electron_bins? "[x]" : "[ ]") <<std::endl
					<<"PE (S2) bins:\t\t" <<(using_S2_bins? "[x]" : "[ ]") <<std::endl;
		if(using_S2_bins || using_S2_threshold)
		{
			std::cout 	<<"\tmu_PE:\t\t" <<S2_mu<<std::endl
			 			<<"\tsigma_PE:\t" <<S2_sigma<<std::endl
			 			<<"\tImported trigger efficiencies:\t" <<(Trigger_Efficiency_PE.empty()? "[ ]" : "[x]") <<std::endl
			 			<<"\tImported acc. efficiencies:\t" <<(Acceptance_Efficiency_PE.empty()? "[ ]" : "[x]") <<std::endl;
			if(using_S2_bins)
			{
				std::cout <<"\n\tBin\tBin range [S2]"<<std::endl;
				for(unsigned int bin = 0; bin < number_of_bins; bin++)
				{
					std::cout<<"\t"<<bin+1<<"\t["<<S2_bin_ranges[bin]<<","<<S2_bin_ranges[bin+1]<<")"<<std::endl;
				}
			}
			
		}
		else if(using_electron_bins || using_electron_threshold)
		{
			std::cout<<"\tNe threshold:\t"<<ne_threshold<<std::endl;
			std::cout<<"\tNe max:\t\t"<<ne_max<<std::endl;
		}
	 	std::cout<<"----------------------------------------"<<std::endl<<std::endl;
	}