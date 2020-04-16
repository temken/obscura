#include "Direct_Detection_Ionization.hpp"

#include <cmath>

//Headers from libphys library
#include "Statistics.hpp"

//1. Event spectra and rates
	// 	double dRdlogEe(double Ee,const DM_Particle& DM,const Atomic_Shell& shell)
// 	{
// 		double mA = shell.Mass_Nucleus;
// 		double qref = aEM*mElectron;
// 		double vDM = 1e-3; //cancels in v^2 dSigmadq2
// 		double prefactor = 1.0/mA*rhoDM/DM.mass/2.0;

// 		return prefactor*sum;
// 	}

	double dRdEe_Ionization(double Ee, const DM_Particle& DM, DM_Distribution& DM_distr, const Atomic_Electron& shell)
	{
		double prefactor = 1.0/shell.nucleus_mass * DM_distr.DM_density / DM.mass / Ee / 2.0;
		double vDM = 1.0e-3; 	//cancels in the product with dSigma_dq^2 
								//-> THIS NEEDS TO BE UPDATED WHEN IMPLEMENTING VELOCITY DEPENDING CROSS SECTIONS
		
		// Integration over q
		// (a) Using numerical integration function.
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
		int ki = std::floor(log10(k / shell.k_min) / shell.dlogk);

		// double vMax = DM_distr.Maximum_DM_Speed();
		// double qMin = DM.mass * vMax - sqrt(DM.mass * DM.mass * vMax * vMax - 2.0 * DM.mass * shell.binding_energy);
		// double qMax = DM.mass * vMax + sqrt(DM.mass * DM.mass * vMax * vMax - 2.0 * DM.mass * shell.binding_energy);
		// int qi_min = std::floor(log10(qMin / shell.q_min) / shell.dlogq);
		// int qi_max = std::floor(log10(qMax / shell.q_min) / shell.dlogq);
		// if(qi_max > shell.Nq) qi_max = shell.Nq;

		int qi_min = 0;
		int qi_max = shell.Nq;
		for(int qi = qi_min; qi < qi_max; qi++)
		{
			double q = shell.q_Grid[qi];
			double vMin = vMinimal_Electrons(q, shell.binding_energy + Ee, DM.mass);
			integral += shell.dlogq * q * q * DM.dSigma_dq2_Electron(q,vDM) * vDM * vDM * DM_distr.Eta_Function(vMin) * shell.Form_Factor_Tables[ki][qi];//shell.Ionization_Form_Factor(q,Ee);
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

	double R_ne_Ionization(unsigned int n_e, const DM_Particle& DM, DM_Distribution& DM_distr, const Atomic_Electron& shell)
	{
		// 		double qref = aEM*mElectron;
// 		double sum=0.0;
// 		for(int ki=0;ki<shell.nk;ki++)
// 		{
// 			double Ee = pow(qref*exp(shell.logkMin+(ki)*shell.dlogk),2.0)/2.0/mElectron;
// 			sum += 2.0 * shell.dlogk * PDFne(ne,Ee,shell) * dRdlogEe(Ee,DM,shell);
// 		}
// 		return sum;
		return 0.0;
	}

	double R_ne_Ionization(unsigned int n_e, const DM_Particle& DM, DM_Distribution& DM_distr, const Atom& atom)
	{
		return 0.0;
	}

	double R_S2_Ionization(unsigned int S2, const DM_Particle& DM, DM_Distribution& DM_distr, const Atomic_Electron& shell)
	{
// 		double sum=0.0;
// 		for(int ne = 1; ne <= 15; ne++)
// 		{
// 				sum += PDF_Gauss(nPE,mu_PE*ne,sqrt(ne)*sigma_PE) * dRdne(ne,DM);
// 		}
// 		if( !(Trigger_Efficiency_PE.empty()) ) sum *= Trigger_Efficiency_PE[nPE-1];
// 		if( !(Acceptance_Efficiency_PE.empty()) ) sum *= Acceptance_Efficiency_PE[nPE-1];
// 		return sum;
		return 0.0;
	}

	double R_S2_Ionization(unsigned int S2, const DM_Particle& DM, DM_Distribution& DM_distr, const Atom& atom)
	{
		return 0.0;
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
		// 		double N=0;
// 		if(Using_Electron_Bins)
// 		{
// 			for(int ne = ne_threshold ; ne < 16 ; ne++)
// 			{
// 				N += dRdne(ne,DM);
// 			}
// 		}
// 		else if(Using_PE_Bins)
// 		{
// 			for(int nPE = bins.front() ; nPE < bins.back() ; nPE++)
// 			{
// 				N += dRdnPE(nPE,DM);
// 			}
// 		}
// 		return N;
		return 0.0;
	}

	void DM_Detector_Ionization::Use_Electron_Bins(unsigned int ne_thr, unsigned int N_bins)
	{
		// 		Using_Electron_Bins = true;
		// 		ne_threshold = n;
	}

	void DM_Detector_Ionization::Use_S2_Bins(unsigned int ne_thr, double S2mu, double S2sigma, unsigned int N_bins)
	{
		// 		Using_PE_Bins = true;
// 		Using_Electron_Bins = false;
// 		mu_PE = mu;
// 		sigma_PE = sigma;
// 		bins = binsizes;
	}

	// 	void Detector_Ionization::Import_Trigger_Efficiency_PE(std::string filename)
// 	{
// 		Trigger_Efficiency_PE = Read_List(filename);
// 	}
// 	void Detector_Ionization::Import_Acceptance_Efficiency_PE(std::string filename)
// 	{
// 		Acceptance_Efficiency_PE = Read_List(filename);
// 	}

	std::vector<double> DM_Detector_Ionization::DM_Signals_Binned(const DM_Particle& DM, DM_Distribution& DM_distr) 
	{
		return {};
	}

	void DM_Detector_Ionization::Print_Summary(int MPI_rank) const 
	{
		// 		Print_Summary_Base();
// 		std::cout 	<<std::endl<<"Electron scattering experiment."<<std::endl
// 					<<"Target:\t\t\t"	<<target_name <<std::endl
// 					<<"Analysis:\t\t" <<(binned_data ? "Binned Poisson" : "Poisson") <<std::endl
// 					<<"PE spectrum:\t\t" <<(Using_PE_Bins? "[x]" : "[ ]") <<std::endl;
// 		if(Using_PE_Bins)
// 		{
// 			std::cout <<"\tmu_PE:\t" <<mu_PE<<std::endl;
// 			std::cout <<"\tsigma_PE:\t" <<sigma_PE<<std::endl;
// 		}
// 		if(binned_data)
// 		{
// 			std::cout <<"Bins:"<<std::endl;
						
// 			if(Using_Electron_Bins)
// 			{
// 				std::cout<<"n_e bin\tevents\tefficiency"<<std::endl;
// 				for(unsigned int bin = 0 ; bin < binned_events.size() ; bin++)
// 				{
// 					std::cout <<bin+1 <<"\t" <<binned_events[bin]<<"\t" <<bin_efficiencies[bin]<<std::endl;
// 				}
// 			}
// 			else if(Using_PE_Bins)
// 			{
// 				std::cout<<"PE bin\t\tevents\tefficiency"<<std::endl;

// 				for(unsigned int bin = 0 ; bin < binned_events.size() ; bin++)
// 				{
// 					std::cout <<"["<<bins[bin]<<","<<bins[bin+1]<<")\t\t" <<binned_events[bin]<<"\t"<<bin_efficiencies[bin]<<std::endl;
// 				}
// 			}
// 		}
// 	 	std::cout<<"----------------------------------------"<<std::endl<<std::endl;
	}