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

	void DM_Detector::Set_Flat_Efficiency(double eff)
	{
		flat_efficiency = eff;
	}

	double DM_Detector::Log_Likelihood(const DM_Particle& DM, DM_Distribution& DM_distr)
	{
		if(statistical_analysis=="Poisson")
		{
			double s = DM_Signals_Total(DM, DM_distr);
			unsigned long int n = observed_events;
			double b = expected_background;
			// if(b < 1.0e-4 && (n > s)) b = n-s;	// see eq.(29) of [arXiv:1705.07920]
			return Log_Likelihood_Poisson(s, n, b);
		}
		else if(statistical_analysis == "Binned Poisson")
		{
			std::vector<double> s = DM_Signals_Binned(DM, DM_distr);
			std::vector<unsigned long int> n = bin_observed_events;
			std::vector<double> b = bin_expected_background;
			// for(unsigned int i = 0; i < b.size(); i++)
			// {
			// 	if(b[i] < 1.0e-4 && (n[i] > s[i])) b[i] = n[i]-s[i]; // see eq.(29) of [arXiv:1705.07920]
			// }
			return Log_Likelihood_Poisson_Binned(s, n, b);
		}
		else if(statistical_analysis == "Maximum-Gap")
		{
			return log(P_Value_Maximum_Gap(DM, DM_distr));
		}
		else
		{
			std::cerr<<"Error in DM_Detector_Nucleus::Log_Likelihood(): Analysis "<<statistical_analysis <<" not recognized."<<std::endl;
			std::exit(EXIT_FAILURE);
		}
	}

	double DM_Detector::Likelihood(const DM_Particle& DM, DM_Distribution& DM_distr)
	{
		return exp( Log_Likelihood(DM, DM_distr) );
	}

	double DM_Detector::P_Value(const DM_Particle& DM, DM_Distribution& DM_distr)
	{
		double p_value = 1.0;
		if(statistical_analysis=="Poisson")
		{
			// //Compute test-statistic t
			// double log_likelihood_0 = Likelihood_Poisson(0.0, observed_events, expected_background );
			// double log_likelihood = Log_Likelihood(DM, DM_distr);
			// double t = -2.0 * (log_likelihood - log_likelihood_0);
			// //Integrate 1/2 chi-square distribution with 1 degree of freedom (dof) over t' larger than t.
			// double dof = 1.0;
			// p_value = 0.5* (1.0 - CDF_Chi_Square(t, dof));

			double expectation_value = DM_Signals_Total(DM, DM_distr) + expected_background;
			p_value = CDF_Poisson(expectation_value, observed_events);
		}
		else if(statistical_analysis == "Binned Poisson")
		{
			// //Compute test-statistic t
			// double log_likelihood_0 = Likelihood_Poisson_Binned(std::vector<double>(number_of_bins,0.0), bin_observed_events, bin_expected_background );
			// double log_likelihood = Log_Likelihood(DM, DM_distr);
			// double t = -2.0 * (log_likelihood - log_likelihood_0);
			// //Integrate chi-square distribution with 1 degree of freedom (dof) over t' larger than t.
			// double dof = 1.0;
			// p_value = 0.5* (1.0 - CDF_Chi_Square(t, dof));

			std::vector<double> expectation_values = DM_Signals_Binned(DM, DM_distr);
			std::vector<double> p_values(number_of_bins, 0.0);
			for(unsigned int i = 0; i < number_of_bins; i++)
			{
				double expectation_value = expectation_values[i] + bin_expected_background[i];
				p_values[i] = CDF_Poisson(expectation_value, bin_observed_events[i]); 
			}
			p_value =  *std::min_element(p_values.begin(),p_values.end());
		}
		else if(statistical_analysis == "Maximum-Gap")
		{
			p_value = P_Value_Maximum_Gap(DM, DM_distr);
		}
		else
		{
			std::cerr<<"Error in DM_Detector_Nucleus::P_Value(): Analysis "<<statistical_analysis <<" not recognized."<<std::endl;
			std::exit(EXIT_FAILURE);
		}
		return p_value;
	}

	double DM_Detector::Upper_Limit(DM_Particle& DM, DM_Distribution& DM_distr, double certainty)
	{
		double interaction_parameter_original = DM.Get_Interaction_Parameter(targets);
		
		std::function<double(double)> func = [this, &DM, &DM_distr, certainty] (double log10_parameter)
		{
			double parameter = pow(10.0, log10_parameter);
			DM.Set_Interaction_Parameter(parameter, targets);
			double p_value = P_Value(DM, DM_distr);
			return p_value-(1.0-certainty);
		};
		double log10_upper_bound = Find_Root(func, -30.0, 10.0, 1.0e-6);

		DM.Set_Interaction_Parameter(interaction_parameter_original, targets);
		return pow(10.0, log10_upper_bound);
	}

	std::vector<std::vector<double>> DM_Detector::Upper_Limit_Curve(DM_Particle& DM, DM_Distribution& DM_distr, double mMin,double mMax, int points, double certainty)
	{
		double mOriginal = DM.mass;
		double lowest_mass = Minimum_DM_Mass(DM, DM_distr);
		std::vector<std::vector<double>> limit;
		std::vector<double> masses = Log_Space(mMin,mMax,points);

		for(unsigned int i = 0; i < masses.size(); i++)
		{
			if(masses[i] < lowest_mass) continue;
			DM.Set_Mass(masses[i]);
			limit.push_back(std::vector<double>{masses[i], Upper_Limit(DM, DM_distr, certainty)});
		}

		DM.Set_Mass(mOriginal);
		return limit;
	}

	//a) Poisson
	void DM_Detector::Set_Observed_Events(unsigned long int n, double B)
	{
		statistical_analysis = "Poisson";
		observed_events = n;

		expected_background = B;
	}

	void DM_Detector::Set_Expected_Background(double B)
	{
		expected_background = B;
	}

	double DM_Detector::DM_Signals_Total(const DM_Particle& DM, DM_Distribution& DM_distr)
	{
		std::function<double(double)> spectrum = [this, &DM, &DM_distr] (double E)
		{
			return dRdE(E, DM, DM_distr);
		};
		double epsilon = Find_Epsilon(spectrum, energy_threshold, energy_max, 1e-6);
		return exposure*Integrate(spectrum, energy_threshold, energy_max, epsilon);
	}

	//b) Binned Poisson
	void DM_Detector::Use_Energy_Bins(double Emin, double Emax, int bins)
	{
		statistical_analysis = "Binned Poisson";
		energy_threshold = Emin;
		energy_max = Emax;

		number_of_bins = bins;
		bin_energies = Linear_Space(energy_threshold, energy_max, bins+1);
		
		if(bin_observed_events.size() != number_of_bins) bin_observed_events = std::vector<unsigned long int>(bins,0);
		if(bin_expected_background.size() != number_of_bins) bin_expected_background = std::vector<double>(bins,0.0);
		if(bin_efficiencies.size() != number_of_bins) bin_efficiencies = std::vector<double>(bins, 1.0);
	}
	
	void DM_Detector::Set_Observed_Events(std::vector<unsigned long int> Ni)
	{
		if(statistical_analysis != "Binned Poisson" || Ni.size() != number_of_bins)
		{
			std::cerr<<"Error in DM_Detector::Set_Observed_Events(std::vector<unsigned long int>): Length of the observed vector is not equal to the number of bins."<<std::endl;
			std::exit(EXIT_FAILURE);
		}
		else
		{
			bin_observed_events = Ni;
			observed_events = std::accumulate(bin_observed_events.begin(),bin_observed_events.end(),0);
		}		
	}

	void DM_Detector::Set_Bin_Efficiencies(const std::vector<double>& eff)
	{
		if(statistical_analysis != "Binned Poisson" || eff.size() != number_of_bins)
		{
			std::cerr<<"Error in DM_Detector::Set_Bin_Efficiencies(const std::vector<double>&): Length of the efficiency vector is not equal to the number of bins."<<std::endl;
			std::exit(EXIT_FAILURE);
		}
		else bin_efficiencies = eff;
	}

	void DM_Detector::Set_Expected_Background(const std::vector<double>& Bi)
	{
		if(statistical_analysis != "Binned Poisson" || Bi.size() != number_of_bins)
		{
			std::cerr<<"Error in DM_Detector::Set_Expected_Background(const std::vector<double>&): Length of the efficiency vector is not equal to the number of bins."<<std::endl;
			std::exit(EXIT_FAILURE);
		}
		else 
		{
			bin_expected_background = Bi;
			expected_background = std::accumulate(bin_expected_background.begin(),bin_expected_background.end(),0);
		}
	}
	

	std::vector<double> DM_Detector::DM_Signals_Binned(const DM_Particle& DM, DM_Distribution& DM_distr)
	{
		if(statistical_analysis != "Binned Poisson" || number_of_bins == 0)
		{
			std::cerr<<"Error in DM_Detector::DM_Signals_Binned(const DM_Particle&, DM_Distribution&): The analysis is not binned Poisson."<<std::endl;
			std::exit(EXIT_FAILURE);
		}
		else if(bin_energies.empty())
		{
			std::cerr<<"Error in DM_Detector::DM_Signals_Binned(const DM_Particle&, DM_Distribution&): No energy bins defined."<<std::endl;
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
				
				double epsilon = Find_Epsilon(spectrum, bin_energies[i], bin_energies[i+1], 1e-4);
				double mu = exposure * Integrate(spectrum, bin_energies[i], bin_energies[i+1], epsilon);
				mu_i.push_back(bin_efficiencies[i] * mu);
			}
			return mu_i;
		}
	}

	void DM_Detector::Use_Maximum_Gap(std::string filename_energy_data,double dim)
	{
		statistical_analysis = "Maximum-Gap";

		maximum_gap_energy_data.clear();
		maximum_gap_energy_data = Import_List(filename_energy_data, dim);
		observed_events = maximum_gap_energy_data.size();
		
		maximum_gap_energy_data.push_back(energy_threshold);
		maximum_gap_energy_data.push_back(energy_max);
		std::sort(maximum_gap_energy_data.begin(),maximum_gap_energy_data.end());
	}

	double DM_Detector::P_Value_Maximum_Gap(const DM_Particle& DM, DM_Distribution& DM_distr)
	{
		// Interpolation spectrum = Spectrum(DM);
		std::function<double(double)> spectrum = [this, &DM, &DM_distr] (double E)
		{
			return exposure * dRdE(E, DM, DM_distr);
		};

		//Determine and save all gap sizes.
		std::vector<double> x;
		for(unsigned int i = 0;i < (maximum_gap_energy_data.size()-1); i++)
		{
			double E1 = maximum_gap_energy_data[i];
			double E2 = maximum_gap_energy_data[i+1];
			double epsilon = Find_Epsilon(spectrum,E1,E2,1e-3);
			double xGap = Integrate(spectrum, E1,E2,epsilon);
			x.push_back(xGap);
		}
		
		//Maximum gap
		double xMax = *std::max_element(x.begin(),x.end());
		
		//P_Value
		double N = std::accumulate(x.begin(),x.end(),0.0);
		double llh = 1.0 - CDF_Maximum_Gap(xMax,N);
		return llh;
	}

	void DM_Detector::Print_Summary_Base(int MPI_rank) const
	{
		if(MPI_rank == 0)
		{
			std::cout 	<<std::endl
						<<"----------------------------------------"<<std::endl
						<<"Experiment summary:\t"<<name<<std::endl
						<<"Target particles:\t" <<targets <<std::endl
						<<"Exposure [kg year]:\t" <<In_Units(exposure,kg*yr)<<std::endl
						<<"Flat efficiency [%]:\t"<<Round(100.0*flat_efficiency)<<std::endl
						<<"Observed events:\t"<<observed_events<<std::endl
						<<"Expected backgroud:\t" <<expected_background <<std::endl
						<<"Statistical analysis:\t" <<statistical_analysis <<std::endl;
			if(statistical_analysis == "Binned Poisson")
			{
				std::cout <<"\tNumber of bins:\t" <<number_of_bins <<std::endl;
				if(!bin_observed_events.empty() && !bin_observed_events.empty())
				{
					std::cout <<"\tBin\tEfficiency[%]\tObserved events\tExpected background"<<std::endl;
					for(unsigned int i = 0; i < bin_observed_events.size(); i++) std::cout <<"\t"<<i+1<<"\t" <<100.0 * bin_efficiencies[i]<<"\t\t"<<bin_observed_events[i] <<"\t\t"<<bin_expected_background[i] <<std::endl;
				}
			}
			std::cout <<std::endl;			
		}
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
