#include "obscura/Target_Crystal.hpp"

#include <cmath>

#include "libphysica/Natural_Units.hpp"
#include "libphysica/Special_Functions.hpp"
#include "libphysica/Utilities.hpp"

#include "version.hpp"

namespace obscura
{
using namespace libphysica::natural_units;

Crystal::Crystal(std::string target)
: N_E(500), N_q(900), name(target), dE(0.1 * eV), dq(0.02 * aEM * mElectron)
{
	E_max = N_E * dE;
	q_max = N_q * dq;

	double prefactor;
	if(name == "Si")
	{
		energy_gap = 1.11 * eV;
		epsilon	   = 3.6 * eV;
		M_cell	   = 2.0 * 28.08 * mNucleon;
		prefactor  = 2.0 * eV;
		Q_max	   = std::floor((E_max - energy_gap + epsilon) / epsilon);
	}
	else if(name == "Ge")
	{
		energy_gap = 0.67 * eV;
		epsilon	   = 2.9 * eV;
		M_cell	   = 2.0 * 72.6 * mNucleon;
		prefactor  = 1.8 * eV;
		Q_max	   = std::floor((E_max - energy_gap + epsilon) / epsilon);
	}
	else
	{
		std::cerr << "Error in obscura::Crystal::Crystal(): target " << target << " not recognized." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	// Import the form factor
	std::string path			 = PROJECT_DIR "data/Semiconductors/C." + target + "137.dat";
	std::vector<double> aux_list = libphysica::Import_List(path);
	std::vector<std::vector<double>> form_factor_table(N_q, std::vector<double>(N_E, 0.0));
	double wk	   = 2.0 / 137.0;
	unsigned int i = 0;
	for(int Ei = 0; Ei < N_E; Ei++)
		for(int qi = 0; qi < N_q; qi++)
			form_factor_table[qi][Ei] = prefactor * (qi + 1) / dE * wk / 4.0 * aux_list[i++];
	std::vector<double> q_grid = libphysica::Linear_Space(dq, q_max, N_q);
	std::vector<double> E_grid = libphysica::Linear_Space(dE, E_max, N_E);
	form_factor_interpolation  = libphysica::Interpolation_2D(q_grid, E_grid, form_factor_table);
}

bool have_warned = false;
double Crystal::Crystal_Form_Factor(double q, double E)
{
	if(q < dq || q > q_max || E < dE || E > E_max)
	{
		if(!have_warned)
		{
			std::cerr << "Warning in obscura::Crystal::Crystal_Form_Factor(): q or E out of range." << std::endl
					  << "\tq = " << libphysica::Round(q / keV) << " keV\tq_min = " << libphysica::Round(dq / keV) << " keV\tq_max = " << libphysica::Round(q_max / keV) << " keV" << std::endl
					  << "\tE = " << libphysica::Round(E / eV) << " eV\tE_min = " << libphysica::Round(dE / eV) << " eV\tE_max = " << libphysica::Round(E_max / eV) << " eV" << std::endl
					  << "\tReturning 0. (This warning will not be repeated.)" << std::endl;
			have_warned = true;
		}
		return 0.0;
	}
	return form_factor_interpolation(q, E);
}

}	// namespace obscura