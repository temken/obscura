#include "Direct_Detection.hpp"

#include <cmath>
#include <functional>
#include <numeric>
#include <algorithm> //for std::min_element, std::max_element, std::sort

//Headers from libphys library
#include "Natural_Units.hpp"
#include "Numerics.hpp"
#include "Statistics.hpp"
#include "Utilities.hpp"

//1. DM_Detector base class

	double DM_Detector::Likelihood_Maximum_Gap(const DM_Particle& DM, DM_Distribution& DM_distr)
	{
		// Interpolation spectrum = Spectrum(DM);
		std::function<double(double)> spectrum = [this, &DM, &DM_distr] (double E)
		{
			return exposure * dRdE(E, DM, DM_distr);
		};

		//Determine and save all gap sizes.
		std::vector<double> x;
		for(unsigned int i = 0;i < (background_energy_data_sorted.size()-1); i++)
		{
			double E1 = background_energy_data_sorted[i];
			double E2 = background_energy_data_sorted[i+1];
			double epsilon = Find_Epsilon(spectrum,E1,E2,1e-3);
			double xGap = Integrate(spectrum, E1,E2,epsilon);
			x.push_back(xGap);
		}
		
		//Maximum gap
		double xMax = *std::max_element(x.begin(),x.end());
		
		//Likelihood
		double N = std::accumulate(x.begin(),x.end(),0.0);
		double llh = 1.0 - CDF_Maximum_Gap(xMax,N);
		return llh;
	}

	void DM_Detector::Print_Summary_Base() const
	{
		std::cout 	<<std::endl
					<<"----------------------------------------"<<std::endl
					<<"Experiment summary:\t"<<name<<std::endl
					<<"Target particles:\t" <<targets <<std::endl
					<<"Exposure [kg day]:\t" <<In_Units(exposure,kg*day)<<std::endl
					<<"Flat efficiency [%]:\t"<<Round(100.0*flat_efficiency)<<std::endl
					<<"Observed events:\t"<<background_events<<std::endl
					<<"Statistical analysis:\t" <<statistical_analysis <<std::endl;
		if(statistical_analysis == "Binned Poisson")
		{
			std::cout <<"\tNumber of bins:\t" <<number_of_bins <<std::endl;
			if(!binned_background.empty())
			{
				std::cout <<"\tBin\tBackground events"<<std::endl;
				for(unsigned int i = 0; i < binned_background.size(); i++) std::cout <<"\t"<<i<<"\t"<<binned_background[i]<<std::endl;
			}
		}
		std::cout <<std::endl;
	}

	void DM_Detector::Set_Flat_Efficiency(double eff)
	{
		flat_efficiency = eff;
	}

	double DM_Detector::Likelihood(const DM_Particle& DM, DM_Distribution& DM_distr)
	{
		double llh;
		if(statistical_analysis=="Poisson")
		{
		 	double expectation_value = Total_Number_of_Signals(DM, DM_distr);
			llh = CDF_Poisson(expectation_value, background_events);
		}
		else if(statistical_analysis == "Binned Poisson")
		{
			std::vector<double> expectation_values = Binned_Number_of_Signals(DM, DM_distr);
			llh = Likelihood_Poisson_Binned(expectation_values, binned_background);
		}
		else if(statistical_analysis == "Maximum-Gap")
		{
			llh = Likelihood_Maximum_Gap(DM, DM_distr);
		}
		else
		{
			std::cerr<<"Error in DM_Detector_Nucleus::Likelihood(): Analysis "<<statistical_analysis <<" not recognized."<<std::endl;
			std::exit(EXIT_FAILURE);
		}
		return llh;
	}

	double DM_Detector::Upper_Bound(DM_Particle& DM, DM_Distribution& DM_distr, double certainty)
	{
		double interaction_parameter_original = DM.Get_Interaction_Parameter();
		
		std::function<double(double)> func = [this, &DM, &DM_distr, certainty] (double log10_parameter)
		{
			double parameter = pow(10.0, log10_parameter);
			DM.Set_Interaction_Parameter(parameter);
			double llh = Likelihood(DM, DM_distr);
			return llh-(1.0-certainty);
		};
		double log10_upper_bound = Find_Root(func, -30.0, 10.0, 1.0e-5);

		DM.Set_Interaction_Parameter(interaction_parameter_original);
		return pow(10.0, log10_upper_bound);
	}

	std::vector<std::vector<double>> DM_Detector::Limit_Curve(DM_Particle& DM, DM_Distribution& DM_distr, double mMin,double mMax, int points, double certainty)
	{
		double mOriginal = DM.mass;
		double lowest_mass = Minimum_DM_Mass(DM, DM_distr);
		std::vector<std::vector<double>> limit(points,std::vector<double>(2,0.0));
		std::vector<double> masses = Log_Space(mMin,mMax,points);

		for(unsigned int i = 0; i < masses.size(); i++)
		{
			if(masses[i] < lowest_mass) continue;
			DM.Set_Mass(masses[i]);
			limit[i][0] = masses[i];
			limit[i][1] = Upper_Bound(DM, DM_distr, certainty);
		}

		DM.Set_Mass(mOriginal);
		return limit;
	}

	//a) Poisson
	void DM_Detector::Set_Background(unsigned long int n)
	{
		statistical_analysis = "Poisson";
		background_events = n;
	}

	double DM_Detector::Total_Number_of_Signals(const DM_Particle& DM, DM_Distribution& DM_distr)
	{
		std::function<double(double)> spectrum = [this, &DM, &DM_distr] (double E)
		{
			return dRdE(E, DM, DM_distr);
		};
		double epsilon = Find_Epsilon(spectrum, energy_threshold, energy_max, 1e-6);
		return exposure*Integrate(spectrum, energy_threshold, energy_max, epsilon);
	}

	//b) Binned Poisson
	void DM_Detector::Define_Energy_Bins(double Emin, double Emax, int bins)
	{
		statistical_analysis = "Binned Poisson";
		number_of_bins = bins;
		bins_energy = Linear_Space(Emin, Emax, bins+1);
		if(binned_background.empty()) binned_background = std::vector<unsigned long int>(bins,0);
	}

	void DM_Detector::Set_Background(std::vector<unsigned long int> Bi)
	{
		if(number_of_bins > 0 && number_of_bins != Bi.size())
		{
			std::cerr<<"Error in DM_Detector::Set_Background(std::vector<unsigned long int>): Length of the background vector is not equal to the number of bins."<<std::endl;
			std::exit(EXIT_FAILURE);
		}
		else
		{
			statistical_analysis = "Binned Poisson";
			number_of_bins = Bi.size();
			binned_background = Bi;
			background_events = std::accumulate(binned_background.begin(),binned_background.end(),0);
		}
		
	}

	std::vector<double> DM_Detector::Binned_Number_of_Signals(const DM_Particle& DM, DM_Distribution& DM_distr)
	{
		if(statistical_analysis != "Binned Poisson" || number_of_bins == 0)
		{
			std::cerr<<"Error in DM_Detector::Binned_Number_of_Signals(const DM_Particle&, DM_Distribution&): The analysis is not binned Poisson."<<std::endl;
			std::exit(EXIT_FAILURE);
		}
		else if(bins_energy.empty())
		{
			std::cerr<<"Error in DM_Detector::Binned_Number_of_Signals(const DM_Particle&, DM_Distribution&): No energy bins defined."<<std::endl;
			std::exit(EXIT_FAILURE);
		}
		else
		{
			// Interpolation spectrum = Spectrum(DM);
			std::function<double(double)> spectrum = [this, &DM, &DM_distr] (double E)
			{
				return dRdE(E, DM, DM_distr);
			};
			std::vector<double> mu_i;
			for(unsigned int i = 0; i < number_of_bins; i++)
			{
				
				double epsilon = Find_Epsilon(spectrum, bins_energy[i], bins_energy[i+1], 1e-4);
				double mu = Integrate(spectrum, bins_energy[i], bins_energy[i+1], epsilon);
				mu_i.push_back(mu);
			}
			return mu_i;
		}
	}

	void DM_Detector::Use_Maximum_Gap(std::string filename_energy_data,double dim)
	{
		statistical_analysis = "Maximum-Gap";

		background_energy_data_sorted.clear();
		background_energy_data_sorted = Import_List(filename_energy_data, dim);
		background_events = background_energy_data_sorted.size();
		
		background_energy_data_sorted.push_back(energy_threshold);
		background_energy_data_sorted.push_back(energy_max);
		std::sort(background_energy_data_sorted.begin(),background_energy_data_sorted.end());
	}




//2. Functions for statistical analysis
	double CDF_Maximum_Gap(double x,double mu)
	{
		if(x==mu) return 1.0-exp(-mu);
		else
		{
			int m = mu/x;
			double sum=0.0;
			for(int k=0;k<=m;k++) 
			{
				double term = pow(k*x-mu,k) / Factorial(k) * exp(-k*x) * (1.0 + k/(mu-k*x));
				sum+= term ;
				if (fabs(term)<1e-20) break;
			}
			return sum;
		}
	}
