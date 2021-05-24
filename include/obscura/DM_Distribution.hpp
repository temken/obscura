#ifndef __DM_Distribution_hpp_
#define __DM_Distribution_hpp_

#include <string>
#include <vector>

#include "libphysica/Linear_Algebra.hpp"
#include "libphysica/Numerics.hpp"

namespace obscura
{

//1. Abstract base class for DM distributions that can be used to compute direct detection recoil spectra.
class DM_Distribution
{
  protected:
	std::string name;
	std::vector<double> v_domain;
	double Eta_Function_Base(double vMin);
	void Print_Summary_Base();

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
	virtual double PDF_Speed(double v);
	virtual double CDF_Speed(double v);

	virtual double Differential_DM_Flux(double v, double mDM);
	virtual double Total_DM_Flux(double mDM);

	//Averages
	virtual libphysica::Vector Average_Velocity();
	virtual double Average_Speed(double vMin = -1.0);

	//Eta-function for direct detection
	virtual double Eta_Function(double vMin);

	virtual void Print_Summary(int mpi_rank = 0);
};

//2. Standard halo model (SHM)
class Standard_Halo_Model : public DM_Distribution
{
  protected:
	double v_0, v_esc, v_observer;
	libphysica::Vector vel_observer;

	//Normalization
	double N_esc;
	virtual void Normalize_PDF();

	double PDF_Velocity_SHM(libphysica::Vector vel);
	double PDF_Speed_SHM(double v);
	double CDF_Speed_SHM(double v);
	double Eta_Function_SHM(double vMin);

	void Print_Summary_SHM();

  public:
	//Constructors:
	Standard_Halo_Model();
	Standard_Halo_Model(double rho, double v0, double vobs, double vesc = 1.0);
	Standard_Halo_Model(double rho, double v0, libphysica::Vector& vel_obs, double vesc = 1.0);

	//Set SHM parameters
	virtual void Set_Speed_Dispersion(double v0);
	void Set_Escape_Velocity(double vesc);
	void Set_Observer_Velocity(const libphysica::Vector& vObserver);
	void Set_Observer_Velocity(int day, int month, int year, int hour = 0, int minute = 0);

	libphysica::Vector Get_Observer_Velocity() const;

	//Distribution functions
	virtual double PDF_Velocity(libphysica::Vector vel) override;
	virtual double PDF_Speed(double v) override;
	virtual double CDF_Speed(double v) override;

	//Eta-function for direct detection
	virtual double Eta_Function(double vMin) override;

	virtual void Print_Summary(int mpi_rank = 0) override;
};

//3. Standard halo model++ (SHM++) as proposed by Evans, O'Hare and McCabe [arXiv:1810.11468]
class SHM_Plus_Plus : public Standard_Halo_Model
{
  protected:
	// Gaia sausage parameters (see the SHM++ paper)
	double eta, beta;
	double sigma_r, sigma_theta, sigma_phi;
	void Compute_Sigmas(double beta);

	//Normalization
	double N_esc_S;
	virtual void Normalize_PDF() override;

	double PDF_Velocity_S(libphysica::Vector vel);
	double PDF_Speed_S(double v);
	double CDF_Speed_S(double v);
	double Eta_Function_S(double vMin);

	// Eta function
	void Interpolate_Eta_Function_S(int v_points = 100);
	libphysica::Interpolation eta_interpolation_s;

	void Print_Summary_SHMpp();

  public:
	//Constructors:
	SHM_Plus_Plus(double e = 0.2, double b = 0.9);
	SHM_Plus_Plus(double rho, double v0, double vobs, double vesc, double e = 0.2, double b = 0.9);
	SHM_Plus_Plus(double rho, double v0, libphysica::Vector& vel_obs, double vesc, double e = 0.2, double b = 0.9);

	//Set SHM++ parameters
	virtual void Set_Speed_Dispersion(double v0) override;

	void Set_Eta(double e);
	void Set_Beta(double b);

	//Distribution functions
	virtual double PDF_Velocity(libphysica::Vector vel) override;
	virtual double PDF_Speed(double v) override;
	virtual double CDF_Speed(double v) override;

	//Eta-function for direct detection
	virtual double Eta_Function(double vMin) override;

	virtual void Print_Summary(int mpi_rank = 0) override;
};

}	// namespace obscura

#endif