#ifndef __DM_Distribution_hpp_
#define __DM_Distribution_hpp_

#include <vector>

#include "Linear_Algebra.hpp"

//1. Abstract base class for DM distributions that can be used to compute direct detection recoil spectra.
class DM_Distribution
{
	protected:
		std::string name;
		
		void Print_Summary_Base() const;
		
	public:
		double DM_density;		//Local DM density
		std::vector<double> v_domain;

		//Constructor:
		DM_Distribution();
		DM_Distribution(std::string label, double rhoDM, double vMax = 1.0, double vMin = 0.0);

		//Distribution functions
		virtual double PDF_Velocity(Vector v) const {return 0.0;};
		virtual double PDF_Speed(double v) const {return 0.0;};
		virtual double CDF_Speed(double v) const;

		virtual Vector Average_Velocity() const;
		virtual double Average_Speed() const;

		//Eta-function for direct detection
		virtual double Eta_Function(double vMin) const;

		virtual void Print_Summary() const {Print_Summary_Base();};
};


#endif