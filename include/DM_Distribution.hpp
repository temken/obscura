#ifndef __DM_Distribution_hpp_
#define __DM_Distribution_hpp_

#include <string>
#include <vector>

//Headers from libphysica library
#include "Linear_Algebra.hpp"

namespace obscura
{

//1. Abstract base class for DM distributions that can be used to compute direct detection recoil spectra.
class DM_Distribution
{
  protected:
	std::string name;
	std::vector<double> v_domain;

	void Print_Summary_Base(int MPI_rank = 0);

  public:
	double DM_density;	 //Local DM density
	bool DD_use_eta_function;

	//Constructors:
	DM_Distribution();
	DM_Distribution(std::string label, double rhoDM, double vMin, double vMax);

	double Minimum_DM_Speed() const;
	double Maximum_DM_Speed() const;

	//Distribution functions
	virtual double PDF_Velocity(libphysica::Vector vel) { return 0.0; };
	virtual double PDF_Speed(double v) { return 0.0; };
	virtual double CDF_Speed(double v);

	virtual double Differential_DM_Flux(double v, double mDM);
	virtual double Total_DM_Flux(double mDM);

	//Averages
	virtual libphysica::Vector Average_Velocity();
	virtual double Average_Speed(double vMin = -1.0);

	//Eta-function for direct detection
	virtual double Eta_Function(double vMin);

	virtual void Print_Summary(int MPI_rank = 0) { Print_Summary_Base(MPI_rank); };
};

//2. Standard halo model (SHM)
class Standard_Halo_Model : public DM_Distribution
{
  protected:
	double v_0, v_esc, v_observer;
	libphysica::Vector vel_observer;

	//Normalization
	double N_esc;
	void Normalize_PDF();

  public:
	//Constructors:
	Standard_Halo_Model();
	Standard_Halo_Model(double rho, double v0, double vobs, double vesc = 1.0);
	Standard_Halo_Model(double rho, double v0, libphysica::Vector& vel_obs, double vesc = 1.0);

	//Set SHM parameters
	void Set_Speed_Dispersion(double v0);
	void Set_Escape_Velocity(double vesc);
	void Set_Observer_Velocity(const libphysica::Vector& vObserver);
	void Set_Observer_Velocity(int day, int month, int year, int hour = 0, int minute = 0);

	//Distribution functions
	virtual double PDF_Velocity(libphysica::Vector vel) override;
	virtual double PDF_Speed(double v) override;

	//Eta-function for direct detection
	virtual double Eta_Function(double vMin) override;

	virtual void Print_Summary(int MPI_rank = 0) override;
};

}	// namespace obscura

#endif