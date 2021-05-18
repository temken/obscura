#include "obscura/Target_Crystal.hpp"

#include <cmath>

#include "libphysica/Natural_Units.hpp"
#include "libphysica/Utilities.hpp"

#include "version.hpp"

namespace obscura
{
using namespace libphysica::natural_units;

Crystal::Crystal(std::string target)
: name(target), dE(0.1 * eV), dq(0.02 * aEM * mElectron)
{
	double prefactor;
	if(name == "Si")
	{
		energy_gap = 1.11 * eV;
		epsilon	   = 3.6 * eV;
		M_cell	   = 2.0 * 28.08 * mNucleon;
		prefactor  = 2.0 * eV;
		Q_max	   = std::floor((50 * eV - energy_gap + epsilon) / epsilon);
	}
	else if(name == "Ge")
	{
		energy_gap = 0.67 * eV;
		epsilon	   = 2.9 * eV;
		M_cell	   = 2.0 * 72.6 * mNucleon;
		prefactor  = 1.8 * eV;
		Q_max	   = std::floor((50 * eV - energy_gap + epsilon) / epsilon);
	}
	else
	{
		std::cerr << "Error in obscura::Crystal::Crystal(): target " << target << " not recognized." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	//Import the form factor
	std::string path			 = PROJECT_DIR "data/Semiconductors/C." + target + "137.dat";
	std::vector<double> aux_list = libphysica::Import_List(path);
	std::vector<std::vector<double>> form_factor_table(900, std::vector<double>(500, 0.0));
	double wk	   = 2.0 / 137.0;
	unsigned int i = 0;
	for(int Ei = 0; Ei < 500; Ei++)
		for(int qi = 0; qi < 900; qi++)
			form_factor_table[qi][Ei] = prefactor * (qi + 1) / dE * wk / 4.0 * aux_list[i++];
	std::vector<double> q_grid = libphysica::Linear_Space(dq, 900 * dq, 900);
	std::vector<double> E_grid = libphysica::Linear_Space(dE, 500 * dE, 500);
	form_factor_interpolation  = libphysica::Interpolation_2D(q_grid, E_grid, form_factor_table);
}

double Crystal::Crystal_Form_Factor(double q, double E)
{
	return form_factor_interpolation(q, E);
}

};	 // namespace obscura