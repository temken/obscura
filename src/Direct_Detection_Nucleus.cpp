#include "Direct_Detection_Nucleus.hpp"

#include <cmath>
#include <numeric> //for std::accumulate
#include <algorithm> //for std::min_element, std::max_element, std::sort

//Headers from libphys library
#include "Statistics.hpp"
#include "Utilities.hpp"

//1. Theoretical nuclear recoil spectrum
	double vMinimal_Nucleus(double ER, double mDM, double mNucleus)
	{
		return sqrt(mNucleus*ER/2.0/pow(Reduced_Mass(mDM,mNucleus),2.0));
	}
	double Maximum_Nuclear_Recoil_Energy(double vDM, double mDM, double mNucleus)
	{
		return 2.0*vDM*vDM*pow(Reduced_Mass(mDM,mNucleus),2.0)/mNucleus;
	}

	double dRdER_Nucleus(double ER, const DM_Particle& DM, DM_Distribution& DM_distr, const Isotope& target_isotope)
	{
		double vMin = vMinimal_Nucleus(ER, DM.mass, target_isotope.mass);
		double rhoDM = DM_distr.DM_density * DM.fractional_density;
		double vDM = 1.0e-3; //cancels
		return 1.0/target_isotope.mass * rhoDM/DM.mass * (vDM*vDM*DM.dSigma_dER_Nucleus(ER,target_isotope,vDM)) * DM_distr.Eta_Function(vMin);
	}

	double dRdER_Nucleus(double ER, const DM_Particle& DM, DM_Distribution& DM_distr, const Element& target_element)
	{
		double dRate=0.0;
	 	for(unsigned int i=0;i<target_element.Number_of_Isotopes();i++)
		{
			dRate+= target_element[i].abundance * dRdER_Nucleus(ER, DM, DM_distr, target_element[i]);
		};
		return dRate;
	}

//2. Nuclear recoil direct detection experiment
	//Constructors
	Detector_Nucleus::Detector_Nucleus()
	: Detector(kg*day,"Nuclei"), target_elements({Get_Element(54)}), relative_mass_fractions({1.0}), energy_threshold(1.0*keV), energy_max(100.0*keV), energy_resolution(0.0), using_efficiency_tables(false), statistical_analysis("Poisson")
	{
	}

	Detector_Nucleus::Detector_Nucleus(double expo,std::vector<Element> elements, double thr, double emax,std::vector<double> abund)
	: Detector(expo,"Nuclei"), target_elements(elements), energy_threshold(thr), energy_max(emax), energy_resolution(0.0), using_efficiency_tables(false), statistical_analysis("Poisson")
	{
		double tot = std::accumulate(abund.begin(),abund.end(),0.0);
		if(abund.empty() || tot > 1.0)
		{
			//Compute relative abundance of the target elements by weight.		
			for(unsigned int i = 0 ; i<target_elements.size() ; i++)
			{
				double proportion = (abund.empty()) ? 1.0 : abund[i];
				double mi = proportion * target_elements[i].Average_Nuclear_Mass();
				relative_mass_fractions.push_back(mi);
			}
			double Mtot = std::accumulate(relative_mass_fractions.begin(),relative_mass_fractions.end(),0.0);
			for(unsigned int i = 0 ; i<relative_mass_fractions.size() ; i++) relative_mass_fractions[i] /= Mtot;
		}
		else relative_mass_fractions = abund;
	}
	
	double Detector_Nucleus::Maximum_Energy_Deposit(const DM_Particle& DM, const DM_Distribution& DM_distr) const
	{
		double vDM = DM_distr.v_domain[1];
		double Emax = 0.0;
		for(unsigned int i = 0 ; i<target_elements.size() ; i++)
		{
			for(unsigned int j = 0 ; j<target_elements[i].Number_of_Isotopes(); j++)
			{
				double ERmax = Maximum_Nuclear_Recoil_Energy(vDM, DM.mass, target_elements[i][j].mass);
				if(ERmax>Emax && DM.Sigma_Nucleus(target_elements[i][j], vDM) > 0.0 ) Emax = ERmax;
			}
		}
		return Emax + 6.0 * energy_resolution;
	}

	double Detector_Nucleus::Minimum_DM_Mass(DM_Particle& DM, const DM_Distribution& DM_distr) const
	{
		std::vector<double> aux;
		double vMax = DM_distr.v_domain[1];
		for(unsigned int i = 0 ; i<target_elements.size() ; i++)
		{
			for(unsigned int j = 0 ; j < target_elements[i].Number_of_Isotopes() ; j++)
			{

				double mMin = target_elements[i][j].mass / (sqrt(2.0*target_elements[i][j].mass/(energy_threshold-2.0*energy_resolution))*vMax-1.0);
				if(DM.Sigma_Nucleus(target_elements[i][j],vMax) > 0.0) aux.push_back(mMin);
			}
		}
		return *std::min_element(aux.begin(),aux.end());
	}

	void Detector_Nucleus::Set_Resolution(double res)
	{
		energy_resolution = res;
	}

	void Detector_Nucleus::Use_Maximum_Gap(std::string filename,double dim)
	{
		statistical_analysis = "Maximum-Gap";

		signal_energies_sorted.clear();
		signal_energies_sorted = Import_List(filename,dim);
		observed_signals = signal_energies_sorted.size();
		
		signal_energies_sorted.push_back(energy_threshold);
		signal_energies_sorted.push_back(energy_max);
		std::sort(signal_energies_sorted.begin(),signal_energies_sorted.end());
	}

	void Detector_Nucleus::Import_Efficiency(std::string filename,double dim)
	{
		using_efficiency_tables = true;
		Interpolation eff(filename,dim);
		efficiencies.push_back(eff);
	}

	void Detector_Nucleus::Import_Efficiency(std::vector<std::string> filenames,double dim)
	{
		efficiencies.clear();
		for(unsigned int i =0; i<filenames.size() ; i++) Import_Efficiency(filenames[i],dim);
	}

	double Detector_Nucleus::dRdE(double E, const DM_Particle& DM, DM_Distribution& DM_distr)
	{

		double dR = 0.0;
		if(energy_resolution < 1e-6*eV)
		{
			double eff = 1.0;
			for(unsigned int i = 0 ; i<target_elements.size() ; i++)
			{
				if(using_efficiency_tables)
				{
					if(efficiencies.size() == 1) eff = efficiencies[0](E);
					else if(efficiencies.size() == target_elements.size()) eff = efficiencies[i](E);
				}
				dR += eff * flat_efficiency * relative_mass_fractions[i] * dRdER_Nucleus(E, DM, DM_distr, target_elements[i]);
			}
		}
		else
		{
			//Find minimum and maximum ER contributing to dR/dE(E):
			std::vector<double> aux = {E - 6.0*energy_resolution, energy_threshold - 2.0*energy_resolution};
			double eMin = *std::max_element(aux.begin(),aux.end());
			double eMax = E + 6.0*energy_resolution;
			
			//Convolute theoretical spectrum with Gaussian
			std::function<double(double)> integrand = [this,E,&DM, &DM_distr] (double ER)
			{
				double dRtheory = 0.0;
				for(unsigned int i = 0 ; i<target_elements.size() ; i++)
				{
					double eff = 1.0;
					if(using_efficiency_tables)
					{
						if(efficiencies.size() == 1) eff = efficiencies[0](E);
						else if(efficiencies.size() == target_elements.size()) eff = efficiencies[i](E);
					}
					dRtheory+= eff * flat_efficiency * relative_mass_fractions[i] * dRdER_Nucleus(ER, DM, DM_distr, target_elements[i]);
				}
				return PDF_Gauss(E,ER,energy_resolution)*dRtheory;
			};
			// double epsilon = 1e-6*(eMax-eMin)*integrand(eMin);
			double epsilon = Find_Epsilon(integrand,eMin,eMax,1e-4);
			dR = Integrate(integrand, eMin,eMax,epsilon);
		}
		return dR;
	}

	double Detector_Nucleus::N_Signals(const DM_Particle& DM, DM_Distribution& DM_distr)
	{
		std::function<double(double)> spectrum = [this, &DM, &DM_distr] (double E)
		{
			return dRdE(E, DM, DM_distr);
		};
		double epsilon = Find_Epsilon(spectrum, energy_threshold, energy_max, 1e-6);
		return exposure*Integrate(spectrum, energy_threshold, energy_max, epsilon);
	}

	double Detector_Nucleus::Likelihood_Maximum_Gap(const DM_Particle& DM, DM_Distribution& DM_distr)
	{
		// Interpolation spectrum = Spectrum(DM);
		std::function<double(double)> spectrum = [this, &DM, &DM_distr] (double E)
		{
			return dRdE(E, DM, DM_distr);
		};
		//Determine and save all gap sizes.
			std::vector<double> x;
			for(unsigned int i = 0;i < (signal_energies_sorted.size()-1); i++)
			{
				double E1 = signal_energies_sorted[i];
				double E2 = signal_energies_sorted[i+1];
				double epsilon = Find_Epsilon(spectrum,E1,E2,1e-3);
				double xGap = exposure*Integrate(spectrum, E1,E2,epsilon);
				x.push_back(xGap);
			}
		//Maximum gap
			double xMax = *std::max_element(x.begin(),x.end());
		//Likelihood
			double N = std::accumulate(x.begin(),x.end(),0.0);
			double llh = 1.0 - CDF_Maximum_Gap(xMax,N);
			return llh;
	}

	double Detector_Nucleus::Likelihood(const DM_Particle& DM, DM_Distribution& DM_distr)
	{
		double llh;
		if(statistical_analysis=="Poisson")
		{
		 	double N = N_Signals(DM, DM_distr);
			llh = CDF_Poisson(N,observed_signals);
		}
		else if(statistical_analysis == "Maximum-Gap")
		{
			llh = Likelihood_Maximum_Gap(DM, DM_distr);
		}
		else
		{
			std::cerr<<"Error in Detector_Nucleus::Likelihood(): Analysis "<<statistical_analysis <<" not recognized."<<std::endl;
			std::exit(EXIT_FAILURE);
		}
		return llh;
	}

	double Detector_Nucleus::Upper_Bound(DM_Particle& DM, DM_Distribution& DM_distr, double certainty)
	{
		double sigma;
		if(statistical_analysis=="Poisson")
		{
		 	double N_Limit = Inv_CDF_Poisson(observed_signals,1.0-certainty);
			double N = N_Signals(DM, DM_distr);
			sigma = N_Limit / N * DM.Sigma_Proton();
		}
		else if(statistical_analysis=="Maximum-Gap")
		{
			double sigma_original = DM.Sigma_Proton();
			double s1=1e-50*cm*cm;
			double s2=1e-25*cm*cm;
			double s3 = pow(10.0, log10(s1*s2)/2.0);
			double e = 1e-4;
			double llh=0.0;
			while(fabs(llh-(1.0-certainty))>e)
			{
				DM.Set_Sigma_Proton(s3);
				llh = Likelihood(DM, DM_distr);
				if(llh>(1.0-certainty))
				{
					s1 = s3;
				}
				else
				{
					s2 = s3;
				}
				s3 = pow(10.0, log10(s1*s2)/2.0);
			}
			sigma=s3;
			DM.Set_Sigma_Proton(sigma_original);
		}
		else
		{
			std::cerr<<"Error in Detector_Nucleus::Upper_Bound(): Analysis "<<statistical_analysis <<" not recognized."<<std::endl;
			std::exit(EXIT_FAILURE);
		}		
		return sigma;
	}

	double Detector_Nucleus::Minimum_DM_Speed(const DM_Particle& DM) const
	{
		double Emin = energy_threshold - 2.0*energy_resolution;
		double vcut = 1.0;
		for(unsigned int i = 0 ; i < target_elements.size() ; i++)
		{
			for(unsigned int j = 0 ; j < target_elements[i].Number_of_Isotopes() ; j++)
			{
				double vmin = vMinimal_Nucleus(Emin, DM.mass, target_elements[i][j].mass);
				if(vmin < vcut && DM.Sigma_Nucleus(target_elements[i][j], 1.0e-3) > 0.0) vcut = vmin;
			}
		}
		return vcut;
	}

	void Detector_Nucleus::Print_Summary() const
	{
		Print_Summary_Base();
		std::cout 	<<std::endl<<"Nuclear recoil experiment." <<std::endl
					<<"Nuclear targets:"	<<std::endl;
		for(unsigned int i=0 ; i < target_elements.size() ; i++)
		{
			std::cout <<"\t" <<target_elements[i].name<<"\t"<<Round(100.0*relative_mass_fractions[i])<<"%"<<std::endl;
			// target_elements[i].Print_Summary();
		}
		std::cout 	<<"Threshold [keV]:\t"<<In_Units(energy_threshold,keV)<<std::endl
					<<"ER_max [keV]:\t\t"<<In_Units(energy_max,keV)<<std::endl
					<<"Energy resolution [eV]:\t"<<In_Units(energy_resolution,eV)<<std::endl
					<<"Analysis:\t\t" <<statistical_analysis <<std::endl
		 			<<"----------------------------------------"<<std::endl<<std::endl;
	}
