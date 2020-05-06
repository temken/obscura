#include "Direct_Detection.hpp"

#include <iostream>
#include <algorithm>
#include <numeric>
#include <cmath>

//Headers from libphysica library
#include "Numerics.hpp"
#include "Statistics.hpp"
#include "Utilities.hpp"

// DM Detector base class, which provides the statistical methods and energy bins.
	//Statistics
	//Likelihoods
	double DM_Detector::Log_Likelihood(const DM_Particle& DM, DM_Distribution& DM_distr)
	{
		if(statistical_analysis == "Poisson")
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
		else if(statistical_analysis == "Maximum Gap")
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
			// p_value = 0.5 * (1.0 - CDF_Chi_Square(t, dof));
			std::vector<double> expectation_values = DM_Signals_Binned(DM, DM_distr);
			std::vector<double> p_values(number_of_bins, 0.0);
			for(unsigned int i = 0; i < number_of_bins; i++)
			{
				double expectation_value = expectation_values[i] + bin_expected_background[i];
				p_values[i] = CDF_Poisson(expectation_value, bin_observed_events[i]); 
			}
			p_value =  *std::min_element(p_values.begin(),p_values.end());
		}
		else if(statistical_analysis == "Maximum Gap")
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

	// (a) Poisson statistics
	void DM_Detector::Initialize_Poisson()
	{
			statistical_analysis = "Poisson";
			observed_events = 0;
			expected_background = 0;
	}

	void DM_Detector::Set_Observed_Events(unsigned long int n)
	{
		if(statistical_analysis != "Poisson")
		{
			std::cerr <<"Error in DM_Detector::Set_Observed_Events(unsigned long int): Statistical analysis is " <<statistical_analysis <<" not 'Poisson'." <<std::endl;
			std::exit(EXIT_FAILURE);
		}
		else
		{
			observed_events = n;
		}
	}

	void DM_Detector::Set_Expected_Background(double B)
	{
		if(statistical_analysis != "Poisson")
		{
			std::cerr <<"Error in DM_Detector::Set_Expected_Background(double): Statistical analysis is " <<statistical_analysis <<" not 'Poisson'." <<std::endl;
			std::exit(EXIT_FAILURE);
		}
		else
		{
			expected_background = B;
		}
	}

	// (b) Binned Poisson statistics
	void DM_Detector::Initialize_Binned_Poisson(unsigned bins)
	{
		if(statistical_analysis == "Binned Poisson")
		{
			std::cerr <<"Error in DM_Detector::Initialize_Binned_Poisson(): Bins have already been defined."<<std::endl;
			std::exit(EXIT_FAILURE);
		}
		else
		{
			statistical_analysis = "Binned Poisson";
			number_of_bins = bins;
			bin_observed_events = std::vector<unsigned long int>(number_of_bins,0);
			bin_expected_background = std::vector<double>(number_of_bins,0.0);
			bin_efficiencies = std::vector<double>(number_of_bins, 1.0);
		}
	}

	void DM_Detector::Set_Observed_Events(std::vector<unsigned long int> Ni)
	{
		if (statistical_analysis != "Binned Poisson")
		{
			std::cerr<<"Error in DM_Detector::Set_Observed_Events(std::vector<unsigned long int>): Statistical analysis is " <<statistical_analysis <<" not 'Binned Poisson'." <<std::endl;
			std::exit(EXIT_FAILURE);
		}
		else if(Ni.size() != number_of_bins)
		{
			std::cerr<<"Error in DM_Detector::Set_Observed_Events(std::vector<unsigned long int>): Length of the input ("<<Ni.size()<<") is not equal to the number of bins(" <<number_of_bins <<")."<<std::endl;
			std::exit(EXIT_FAILURE);
		}
		else
		{
			bin_observed_events = Ni;
			observed_events = std::accumulate(bin_observed_events.begin(), bin_observed_events.end(), 0);
		}		
	}

	void DM_Detector::Set_Bin_Efficiencies(const std::vector<double>& eff)
	{
		if (statistical_analysis != "Binned Poisson")
		{
			std::cerr<<"Error in DM_Detector::Set_Bin_Efficiencies(const std::vector<double>&): Statistical analysis is " <<statistical_analysis <<" not 'Binned Poisson'." <<std::endl;
			std::exit(EXIT_FAILURE);
		}
		else if(eff.size() != number_of_bins)
		{
			std::cerr<<"Error in DM_Detector::Set_Bin_Efficiencies(const std::vector<double>&): Length of the input ("<<eff.size()<<") is not equal to the number of bins(" <<number_of_bins <<")."<<std::endl;
			std::exit(EXIT_FAILURE);
		}
		else
		{
			bin_efficiencies = eff;
		}	
	}

	void DM_Detector::Set_Expected_Background(const std::vector<double>& Bi)
	{
		if (statistical_analysis != "Binned Poisson")
		{
			std::cerr<<"Error in DM_Detector::Set_Expected_Background(const std::vector<double>&): Statistical analysis is " <<statistical_analysis <<" not 'Binned Poisson'." <<std::endl;
			std::exit(EXIT_FAILURE);
		}
		else if(Bi.size() != number_of_bins)
		{
			std::cerr<<"Error in DM_Detector::Set_Expected_Background(const std::vector<double>&): Length of the input ("<<Bi.size()<<") is not equal to the number of bins(" <<number_of_bins <<")."<<std::endl;
			std::exit(EXIT_FAILURE);
		}
		else
		{
			bin_expected_background = Bi;
			expected_background = std::accumulate(bin_expected_background.begin(), bin_expected_background.end(), 0.0);
		}	
	}
	
	// (c) Maximum gap a'la Yellin
	void DM_Detector::Use_Maximum_Gap(std::vector<double> energies)
	{
		statistical_analysis = "Maximum Gap";

		maximum_gap_energy_data = energies;
		std::sort(maximum_gap_energy_data.begin(),maximum_gap_energy_data.end());
		energy_threshold = maximum_gap_energy_data.front();
		energy_max = maximum_gap_energy_data.back();
	}

	double CDF_Maximum_Gap(double x,double mu)
	{
		if(x == mu) return 1.0 - exp(-mu);
		else
		{
			int m = mu/x;
			double sum=0.0;
			for(int k=0; k <= m; k++) 
			{
				double term = pow(k*x-mu,k) / Factorial(k) * exp(-k*x) * (1.0 + k/(mu-k*x));
				sum += term;
				if(fabs(term) < 1e-20) break;
			}
			return sum;
		}
	}

	double DM_Detector::P_Value_Maximum_Gap(const DM_Particle& DM, DM_Distribution& DM_distr)
	{
		std::function<double(double)> spectrum = [this, &DM, &DM_distr] (double E)
		{
			return exposure * dRdE(E, DM, DM_distr);
		};

		//Determine all gaps and find the maximum.
		std::vector<double> gaps;
		for(unsigned int i = 0; i < (maximum_gap_energy_data.size()-1); i++)
		{
			double E1 = maximum_gap_energy_data[i];
			double E2 = maximum_gap_energy_data[i+1];
			double eps = Find_Epsilon(spectrum,E1,E2,1e-3);
			double gap = Integrate(spectrum, E1,E2,eps);
			gaps.push_back(gap);
		}
		
		double max_gap = *std::max_element(gaps.begin(),gaps.end());
		
		double N = std::accumulate(gaps.begin(),gaps.end(),0.0);
		double p_value = 1.0 - CDF_Maximum_Gap(max_gap,N);
		return p_value;
	}


	void DM_Detector::Set_Flat_Efficiency(double eff)
	{
		flat_efficiency = eff;
	}

	//DM functions
	double DM_Detector::DM_Signals_Total(const DM_Particle& DM, DM_Distribution& DM_distr)
	{
		double N=0;
		if(statistical_analysis == "Binned Poisson")
		{
			std::vector<double> binned_events = DM_Signals_Binned(DM, DM_distr);
			N = std::accumulate(binned_events.begin(), binned_events.end(), 0.0);
		}
		else
		{
			std::function<double(double)> spectrum = [this, &DM, &DM_distr] (double E)
			{
				return dRdE(E, DM, DM_distr);
			};
			double epsilon = Find_Epsilon(spectrum, energy_threshold, energy_max, 1e-6);
			N = exposure * Integrate(spectrum, energy_threshold, energy_max, epsilon);
		}
		return N;
	}

	std::vector<double> DM_Detector::DM_Signals_Binned(const DM_Particle& DM, DM_Distribution& DM_distr)
	{
		if(statistical_analysis != "Binned Poisson")
		{
			std::cerr<<"Error in DM_Detector::DM_Signals_Binned(const DM_Particle&, DM_Distribution&): The statistical analysis is " <<statistical_analysis <<", not 'Binned Poisson'."<<std::endl;
			std::exit(EXIT_FAILURE);
		}
		else if(using_energy_bins)
		{
			return DM_Signals_Energy_Bins(DM, DM_distr);
		}
		else
		{
			std::cerr <<"Error in DM_Detector::DM_Signals_Binned(): Statistical analysis is 'Binned Poisson' but no bins have been defined. This should not happen ever." <<std::endl;
			std::exit(EXIT_FAILURE);
		}
	}

	//Limits/Constraints
	double DM_Detector::Upper_Limit(DM_Particle& DM, DM_Distribution& DM_distr, double certainty)
	{
		double interaction_parameter_original = DM.Get_Interaction_Parameter(targets);
		// Find the interaction parameter such that p = 1-certainty
		std::function<double(double)> func = [this, &DM, &DM_distr, certainty] (double log10_parameter)
		{
			double parameter = pow(10.0, log10_parameter);
			DM.Set_Interaction_Parameter(parameter, targets);
			double p_value = P_Value(DM, DM_distr);
			return p_value - (1.0-certainty);
		};
		double log10_upper_bound = Find_Root(func, -30.0, 10.0, 1.0e-4);

		DM.Set_Interaction_Parameter(interaction_parameter_original, targets);
		return pow(10.0, log10_upper_bound);
	}

	std::vector<std::vector<double>> DM_Detector::Upper_Limit_Curve(DM_Particle& DM, DM_Distribution& DM_distr, std::vector<double> masses, double certainty)
	{
		double mOriginal = DM.mass;
		double lowest_mass = Minimum_DM_Mass(DM, DM_distr);
		std::vector<std::vector<double>> limit;

		for(unsigned int i = 0; i < masses.size(); i++)
		{
			if(masses[i] < lowest_mass) continue;
			DM.Set_Mass(masses[i]);
			limit.push_back(std::vector<double>{masses[i], Upper_Limit(DM, DM_distr, certainty)});
			std::cout 	<<i+1 <<"/"<<masses.size()
						<<"\tmDM = "<<Round(In_Units(DM.mass, (DM.mass<GeV)? MeV : GeV)) <<((DM.mass<GeV)? " MeV" : " GeV")
						<<"\tUpper Bound:\t" <<Round(In_Units(limit.back()[1],cm*cm)) <<std::endl;
		}
		DM.Set_Mass(mOriginal);
		return limit;
	}

	//Energy spectrum
	void DM_Detector::Use_Energy_Threshold(double Ethr, double Emax)
	{
		Initialize_Poisson();
		using_energy_threshold = true;
		energy_threshold = Ethr;
		energy_max = Emax;
		if(energy_max < energy_threshold)
		{
			std::cerr <<"Error in DM_Detector::Use_Energy_Threshold(): Energy threshold (" <<energy_threshold/keV <<"keV) is higher than maximum energy (" <<energy_max/keV<<"keV)."<<std::endl;
			std::exit(EXIT_FAILURE);
		}
	}

	void DM_Detector::Use_Energy_Bins(double Emin, double Emax, int bins)
	{
		Initialize_Binned_Poisson(bins);
		using_energy_bins = true;
		energy_threshold = Emin;
		energy_max = Emax;
		bin_energies = Linear_Space(energy_threshold, energy_max, number_of_bins + 1);
		if(energy_max < energy_threshold)
		{
			std::cerr <<"Error in DM_Detector::Use_Energy_Bins(): Energy threshold (" <<energy_threshold/keV <<"keV) is higher than maximum energy (" <<energy_max/keV<<"keV)."<<std::endl;
			std::exit(EXIT_FAILURE);
		}
	}

	std::vector<double> DM_Detector::DM_Signals_Energy_Bins(const DM_Particle& DM, DM_Distribution& DM_distr)
	{
		if(!using_energy_bins)
		{
			std::cerr <<"Error in DM_Detector::DM_Signals_Energy_Bins(const DM_Particle&,DM_Distribution&): Not using energy bins." <<std::endl;
			std::exit(EXIT_FAILURE);
		}
		else
		{
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


	void DM_Detector::Print_Summary_Base(int MPI_rank) const
	{
		if(MPI_rank == 0)
		{
			std::cout 	<<std::endl
						<<"----------------------------------------"<<std::endl
						<<"Experiment summary:\t"<<name<<std::endl
						<<"\tTarget particles:\t" <<targets <<std::endl
						<<"\tExposure [kg year]:\t" <<Round(In_Units(exposure,kg*yr))<<std::endl
						<<"\tFlat efficiency [%]:\t"<<Round(100.0*flat_efficiency)<<std::endl
						<<"\tObserved events:\t"<<observed_events<<std::endl
						<<"\tExpected background:\t" <<expected_background <<std::endl
						<<"\tStatistical analysis:\t" <<statistical_analysis <<std::endl;
			if(statistical_analysis == "Binned Poisson")
			{
				std::cout <<"\t\tNumber of bins:\t" <<number_of_bins <<std::endl;
				if(!bin_observed_events.empty() && !bin_observed_events.empty())
				{
					std::cout <<"\t\tBin\tEfficiency[%]\tObserved events\tExpected background"<<std::endl;
					for(unsigned int i = 0; i < bin_observed_events.size(); i++) std::cout <<"\t\t"<<i+1<<"\t" <<100.0 * bin_efficiencies[i]<<"\t\t"<<bin_observed_events[i] <<"\t\t"<<bin_expected_background[i] <<std::endl;
				}
			}
			if(using_energy_threshold || using_energy_bins || statistical_analysis == "Maximum Gap")
				std::cout <<"\tRecoil energies [keV]:\t["<<Round(energy_threshold/keV)<<","<<Round(energy_max/keV) <<"]"<<std::endl;
			if(using_energy_bins)
			{
				std::cout <<"\n\t\tBin\tBin range [keV]"<<std::endl;
				for(unsigned int bin = 0; bin < number_of_bins; bin++)
				{
					std::cout<<"\t\t"<<bin+1<<"\t["<<Round(In_Units(bin_energies[bin],keV))<<","<<Round(In_Units(bin_energies[bin+1],keV))<<")"<<std::endl;
				}
			}
			std::cout <<std::endl;			
		}
	}
