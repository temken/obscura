#include "obscura/Target_Crystal.hpp"

#include <cmath>

#include "libphysica/Natural_Units.hpp"
#include "libphysica/Special_Functions.hpp"
#include "libphysica/Statistics.hpp"
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

		// Import the ionization yield tables
		auto raw_data = libphysica::Import_Table(PROJECT_DIR "data/Semiconductors/Ionization_yield_Silicon_p100K.dat", {eV, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0});
		for(int Q = 1; Q <= 20; Q++)
		{
			std::vector<double> E_grid, yield_grid;
			for(int i = 0; i < raw_data.size(); i++)
			{
				E_grid.push_back(raw_data[i][0]);
				yield_grid.push_back(raw_data[i][Q]);
			}
			ionization_yield_interpolations.push_back(libphysica::Interpolation(E_grid, yield_grid));
		}
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
		std::cerr << libphysica::Formatted_String("Error", "Red", true) << " in obscura::Crystal::Crystal(): target " << target << " not recognized." << std::endl;
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

double Crystal::Ionization_Yield(double Ee, unsigned int Q)
{
	if(Ee < energy_gap)
		return 0.0;
	else if(name == "Si")
	{
		if(Ee >= ionization_yield_interpolations[Q - 1].domain[0] && Ee <= ionization_yield_interpolations[Q - 1].domain[1])
			return ionization_yield_interpolations[Q - 1](Ee);
		else
		{
			// Above 50 eV use a Gaussian. See Eq. (15) of https://journals.aps.org/prd/pdf/10.1103/PhysRevD.102.063026
			double epsilon_eh_inf = epsilon;   // mean energy per electron-hole pair
			double F_inf		  = 0.119;	   // Fano factor from 2004.11499

			double mean	 = Q * epsilon_eh_inf;
			double sigma = std::sqrt(Q * F_inf) * epsilon_eh_inf;
			return epsilon_eh_inf * libphysica::PDF_Gauss(Ee, mean, sigma);
		}
	}
	else if(name == "Ge")
	{
		double E_1 = epsilon * (Q - 1) + energy_gap;
		double E_2 = epsilon * Q + energy_gap;
		if(Ee < E_1 || Ee > E_2)
			return 0.0;
		else
			return 1.0;
	}
	else
	{
		libphysica::Check_For_Error(true, "obscura::Crystal::Ionization_Yield()", "Target " + name + " not recognized.");
		std::exit(EXIT_FAILURE);
	}
}

bool have_warned = false;
double Crystal::Crystal_Form_Factor(double q, double E)
{
	if(q < dq || q > q_max || E < dE || E > E_max)
	{
		if(!have_warned && (q < 0.999999 * dq || q > 1.000001 * q_max || E < 0.999999 * dE || E > 1.000001 * E_max))
		{
			std::cerr << libphysica::Formatted_String("Warning", "Yellow", true) << " in obscura::Crystal::Crystal_Form_Factor(): q or E out of range." << std::endl
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