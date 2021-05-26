#ifndef __DM_Distribution_hpp_
#define __DM_Distribution_hpp_

#include <string>
#include <vector>

#include "libphysica/Linear_Algebra.hpp"
#include "libphysica/Numerics.hpp"

namespace obscura
{

// 1. Abstract base class for DM distributions that can be used to compute direct detection recoil spectra.
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
	virtual double PDF_Norm();

	virtual double Differential_DM_Flux(double v, double mDM);
	virtual double Total_DM_Flux(double mDM);

	//Averages
	virtual libphysica::Vector Average_Velocity();
	virtual double Average_Speed(double vMin = -1.0);

	//Eta-function for direct detection
	virtual double Eta_Function(double vMin);

	virtual void Print_Summary(int mpi_rank = 0);
	void Export_PDF_Speed(std::string file_path, int v_points = 100, bool log_scale = false);
	void Export_Eta_Function(std::string file_path, int v_points = 100, bool log_scale = false);
};

// 2. Import a tabulated DM distribution from a file (format v[km/sec] :: f(v) [sec/km])
class Imported_DM_Distribution : public DM_Distribution
{
  protected:
	std::string file_path;
	libphysica::Interpolation pdf_speed, eta_function;

	void Check_Normalization();

	double Eta_Function_Int(double v_min);
	void Interpolate_Eta();

  public:
	Imported_DM_Distribution(double rho, const std::string& filepath);

	virtual double PDF_Speed(double v) override;
	virtual double CDF_Speed(double v) override;

	virtual double Eta_Function(double vMin) override;

	virtual void Print_Summary(int mpi_rank = 0) override;
};
}	// namespace obscura

#endif