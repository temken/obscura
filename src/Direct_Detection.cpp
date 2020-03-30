#include "Direct_Detection.hpp"

//Headers from libphys library
#include "Natural_Units.hpp"
#include "Utilities.hpp"

//1. Detector base class
	// Interpolation Detector::Spectrum_Base(const DM_Particle& DM, const DM_Distribution& DM_distr, double Emin,double Emax,unsigned int points)
	// {
	// 	//1. Find maximum ER, i.e. the domain of the spectrum
	// 		double ERmax = Maximum_Energy_Deposit(DM);
	// 		ERmax = std::min(ERmax,Emax);
	// 	//2. Tabulate the spectrum
	// 		double dER = (ERmax-Emin)/(points-1.0);
	// 		std::vector<std::vector<double>> interpol_list;
	// 		for(unsigned int i = 0;i<points;i++)
	// 		{
	// 			double ER = Emin + i*dER;
	// 			double dR = dRdE(ER,DM);
	// 			interpol_list.push_back(std::vector<double> {ER,dR});
	// 		}
	// 	//3. if the maximum value is not Emax we append a last point
	// 		if(ERmax<Emax) interpol_list.push_back(std::vector<double> {Emax,0.0});
	// 	//4. Interpolate and return
	// 		Interpolation spectrum(interpol_list);
	// 		return spectrum;
	// }

	void Detector::Print_Summary_Base() const
	{
		std::cout 	<<std::endl
					<<"----------------------------------------"<<std::endl
					<<"Experiment summary:\t"<<name<<std::endl
					<<"Target particles:\t" <<targets <<std::endl
					<<"Exposure [kg day]:\t" <<In_Units(exposure,kg*day)<<std::endl
					<<"Flat efficiency [%]:\t"<<Round(100.0*flat_efficiency)<<std::endl
					<<"Observed events:\t"<<observed_signals<<std::endl;
	}

	void Detector::Set_Name(std::string n)
	{
		name=n;
	}

	void Detector::Set_Flat_Efficiency(double eff)
	{
		flat_efficiency = eff;
	}

	void Detector::Set_Observed_Signals(unsigned long int n)
	{
		observed_signals = n;
	}

	std::vector<std::vector<double>> Detector::Limit_Curve(DM_Particle& DM, const DM_Distribution& DM_distr, double mMin,double mMax, int points, double certainty)
	{
		double mOriginal = DM.mass;
		std::vector<std::vector<double>> limit(points,std::vector<double>(2,0.0));
		std::vector<double> masses = Log_Space(mMin,mMax,points);

		for(unsigned int i = 0; i < masses.size(); i++)
		{
			DM.Set_Mass(masses[i]);
			limit[i][0] = masses[i];
			limit[i][1] = Upper_Bound(DM, DM_distr, certainty);
		}

		DM.Set_Mass(mOriginal);
		return limit;
	}
