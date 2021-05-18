#ifndef __Direct_Detection_Nucleus_hpp_
#define __Direct_Detection_Nucleus_hpp_

#include <string>
#include <vector>

#include "libphysica/Natural_Units.hpp"
#include "libphysica/Numerics.hpp"

#include "obscura/DM_Distribution.hpp"
#include "obscura/DM_Particle.hpp"
#include "obscura/Direct_Detection.hpp"

namespace obscura
{

//1. Theoretical nuclear recoil spectrum [events per time, energy, and target mass]
extern double dRdER_Nucleus(double ER, const DM_Particle& DM, DM_Distribution& DM_distr, const Isotope& target_isotope);
extern double dRdER_Nucleus(double ER, const DM_Particle& DM, DM_Distribution& DM_distr, const Nucleus& target_nucleus);

//2. Nuclear recoil direct detection experiment
class DM_Detector_Nucleus : public DM_Detector
{
  private:
	//Target material
	std::vector<Nucleus> target_nuclei;
	std::vector<double> relative_mass_fractions;

	//Experimental parameters
	double energy_resolution;
	bool using_efficiency_tables;
	std::vector<libphysica::Interpolation> efficiencies;

	virtual double Maximum_Energy_Deposit(const DM_Particle& DM, const DM_Distribution& DM_distr) const override;

  public:
	DM_Detector_Nucleus();
	DM_Detector_Nucleus(std::string label, double expo, std::vector<Nucleus> nuclei, std::vector<double> abund = {});

	void Set_Resolution(double res);
	void Import_Efficiency(std::string filename, double dim);
	void Import_Efficiency(std::vector<std::string> filenames, double dim);

	virtual double Minimum_DM_Speed(const DM_Particle& DM) const override;
	virtual double Minimum_DM_Mass(DM_Particle& DM, const DM_Distribution& DM_distr) const override;
	virtual double dRdE(double E, const DM_Particle& DM, DM_Distribution& DM_distr) override;

	virtual void Print_Summary(int MPI_rank = 0) const override;
};

}	// namespace obscura

#endif