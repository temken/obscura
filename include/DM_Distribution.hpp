#ifndef __DM_Distribution_hpp_
#define __DM_Distribution_hpp_

#include <vector>

//Headers from libphys library
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

			//Constructors:
			DM_Distribution();
			DM_Distribution(std::string label, double rhoDM, double vMin, double vMax);

			//Distribution functions
			virtual double PDF_Velocity(Vector vel) const {return 0.0;};
			virtual double PDF_Speed(double v) const {return 0.0;};
			virtual double CDF_Speed(double v) const;

			virtual Vector Average_Velocity() const;
			virtual double Average_Speed() const;

			//Eta-function for direct detection
			virtual double Eta_Function(double vMin) const;

			virtual void Print_Summary() const {Print_Summary_Base();};
	};

//2. Standard halo model (SHM)
	class Standard_Halo_Model : public DM_Distribution
	{
		protected:
			double v_0, v_esc, v_observer;
			Vector vel_observer;
			double N_esc;
		public:
			//Constructors:
			Standard_Halo_Model();
			Standard_Halo_Model(double rho, double v0, double vobs, double vesc = 1.0);
			Standard_Halo_Model(double rho, double v0, Vector& vel_obs, double vesc = 1.0);

			//Set SHM parameters
			void Set_Speed_Dispersion(double v0);
			void Set_Escape_Velocity(double vesc);
			void Set_Observer_Velocity(Vector& vObserver);
			void Set_Observer_Velocity(int day, int month, int year, int hour = 0, int minute = 0);

			//Compute N_esc
			void Normalize_PDF();

			//Distribution functions
			virtual double PDF_Velocity(Vector vel) const override;
			virtual double PDF_Speed(double v) const override;

			//Eta-function for direct detection
			virtual double Eta_Function(double vMin) const override;

			virtual void Print_Summary() const override;
	};


#endif