#include "obscura/Target_Atom.hpp"

#include <cmath>
#include <fstream>

#include "libphysica/Natural_Units.hpp"
#include "libphysica/Utilities.hpp"

namespace obscura
{
using namespace libphysica::natural_units;

// 1. Kinematic functions
double vMinimal_Electrons(double q, double Delta_E, double mDM)
{
	return (Delta_E / q + q / 2.0 / mDM);
}

// 3. Bound electrons in isolated atoms
std::string s_names[5] = {"s", "p", "d", "f", "g"};

Atomic_Electron::Atomic_Electron(std::string element, int N, int L, double Ebinding, double kMin, double kMax, double qMin, double qMax, unsigned int neSecondary)
: k_min(kMin), k_max(kMax), q_min(qMin), q_max(qMax), n(N), l(L), binding_energy(Ebinding), number_of_secondary_electrons(neSecondary)
{
	name = element + "_" + std::to_string(n) + s_names[l];
	// Import the tables.
	for(int response = 1; response <= 1; response++)   // NEEDS TO BE PUT BACK TO 4 AFTER THE RESPONSE TABLES 2-4 ARE ADDED
	{
		std::string path									= PROJECT_DIR "data/Atomic_Response_Functions/" + name + "_" + std::to_string(response) + ".txt";
		std::vector<std::vector<double>> form_factor_tables = libphysica::Import_Table(path, {}, 3);
		Nk													= form_factor_tables.size();
		Nq													= form_factor_tables[0].size();
		k_Grid												= libphysica::Log_Space(k_min, k_max, Nk);
		q_Grid												= libphysica::Log_Space(q_min, q_max, Nq);
		atomic_response_interpolations.push_back(libphysica::Interpolation_2D(k_Grid, q_Grid, form_factor_tables));
		dlogk = log10(k_max / k_min) / (Nk - 1.0);
		dlogq = log10(q_max / q_min) / (Nq - 1.0);
	}
}

double Atomic_Electron::Atomic_Response_Function(int response, double q, double E)
{
	double k = sqrt(2.0 * mElectron * E);
	if(q > 1.000001 * q_max || k > 1.000001 * k_max || k < 0.999999 * k_min)
	{
		if(!have_warned)
		{
			std::cerr << libphysica::Formatted_String("Warning", "Yellow", true) << " in Atomic_Response_Function(): Arguments of response " << response << " of " << name << " are out of bound." << std::endl;
			if(q > 1.000001 * q_max)
				std::cerr << "\tq = " << q / keV << " keV\ttabulated q domain: [" << q_min / keV << ", " << q_max / keV << "] keV" << std::endl;
			if(k > 1.000001 * k_max || k < 0.999999 * k_min)
				std::cerr << "\tk = " << k / keV << " keV\ttabulated k domain: [" << k_min / keV << ", " << k_max / keV << "] keV" << std::endl;
			std::cerr << "\tReturning 0. (This warning will not be repeated for " << name << ".)" << std::endl;
			have_warned = true;
		}
		return 0.0;
	}
	else if(response == 1)
	{
		if(q > q_min)
			return atomic_response_interpolations[0](k, q);
		else
		{
			// Dipole approximation for low q
			// See eq. 6 of arXiv:1908.10881
			double q_0	= q_min;
			double FF_0 = atomic_response_interpolations[0](k, q_0);
			return q * q / q_0 / q_0 * FF_0;
		}
	}
	else if(response == 2 || response == 3 || response == 4)
	{
		if(q > q_min)
			return atomic_response_interpolations[response - 1](k, q);
		else
		{
			if(!have_warned)
			{
				std::cerr << libphysica::Formatted_String("Warning", "Yellow", true) << " in Atomic_Response_Function(): Arguments of response " << response << " of " << name << " are out of bound." << std::endl
						  << "\tq = " << q / keV << " keV\ttabulated q domain: [" << q_min / keV << ", " << q_max / keV << "] keV" << std::endl
						  << "\tReturning 0. (This warning will not be repeated for " << name << ".)" << std::endl;
				have_warned = true;
			}
			return 0.0;
		}
	}
	else
	{
		std::cerr << libphysica::Formatted_String("Error", "Red", true) << " in Atomic_Response_Function(): Invalid response function " << response << " of " << name << "." << std::endl;
		std::exit(EXIT_FAILURE);
	}
}

double Atomic_Electron::Ionization_Form_Factor(double q, double E)
{
	return Atomic_Response_Function(1, q, E);
}

void Atomic_Electron::Print_Summary(unsigned int MPI_rank) const
{
	if(MPI_rank == 0)
		std::cout << name << "\t" << In_Units(binding_energy, eV) << "\t\t" << number_of_secondary_electrons << std::endl;
}

Atom::Atom(int z, double w, std::vector<Atomic_Electron> shells)
: nucleus(Get_Nucleus(z)), electrons(shells), W(w)
{
}

Atom::Atom(const std::string& element_name)
: nucleus(Get_Nucleus(element_name))
{
	if(element_name == "Xe" || element_name == "Xenon")
	{
		double q_min = 1.0 * keV;
		double q_max = 1000.0 * keV;
		double k_min = 0.1 * keV;
		double k_max = 500.0 * keV;
		Atomic_Electron Xe_1s("Xe", 1, 0, 33317.6 * eV, k_min, k_max, q_min, q_max, 0);
		Atomic_Electron Xe_2s("Xe", 2, 0, 5149.21 * eV, k_min, k_max, q_min, q_max, 0);
		Atomic_Electron Xe_2p("Xe", 2, 1, 4837.71 * eV, k_min, k_max, q_min, q_max, 0);
		Atomic_Electron Xe_3s("Xe", 3, 0, 1093.24 * eV, k_min, k_max, q_min, q_max, 0);
		Atomic_Electron Xe_3p("Xe", 3, 1, 958.43 * eV, k_min, k_max, q_min, q_max, 0);
		Atomic_Electron Xe_3d("Xe", 3, 2, 710.73 * eV, k_min, k_max, q_min, q_max, 0);
		Atomic_Electron Xe_4s("Xe", 4, 0, 213.781 * eV, k_min, k_max, q_min, q_max, 3);
		Atomic_Electron Xe_4p("Xe", 4, 1, 163.495 * eV, k_min, k_max, q_min, q_max, 6);
		Atomic_Electron Xe_4d("Xe", 4, 2, 75.5897 * eV, k_min, k_max, q_min, q_max, 4);
		Atomic_Electron Xe_5s("Xe", 5, 0, 25.6986 * eV, k_min, k_max, q_min, q_max, 0);
		Atomic_Electron Xe_5p("Xe", 5, 1, 12.4433 * eV, k_min, k_max, q_min, q_max, 0);

		W		  = 13.8 * eV;
		electrons = {Xe_5p, Xe_5s, Xe_4d, Xe_4p, Xe_4s, Xe_3d, Xe_3p, Xe_3s, Xe_2p, Xe_2s, Xe_1s};
	}
	else if(element_name == "Ar" || element_name == "Argon")
	{
		double q_min = 1.0 * keV;
		double q_max = 1000.0 * keV;
		double k_min = 0.1 * keV;
		double k_max = 100.0 * keV;
		Atomic_Electron Ar_1s("Ar", 1, 0, 3227.55 * eV, k_min, k_max, q_min, q_max, 0);
		Atomic_Electron Ar_2s("Ar", 2, 0, 335.303 * eV, k_min, k_max, q_min, q_max, 0);
		Atomic_Electron Ar_2p("Ar", 2, 1, 260.453 * eV, k_min, k_max, q_min, q_max, 0);
		Atomic_Electron Ar_3s("Ar", 3, 0, 34.7585 * eV, k_min, k_max, q_min, q_max, 0);
		Atomic_Electron Ar_3p("Ar", 3, 1, 16.0824 * eV, k_min, k_max, q_min, q_max, 0);

		W		  = 19.6 * eV;
		electrons = {Ar_3p, Ar_3s, Ar_2p, Ar_2s, Ar_1s};
	}
	else
	{
		std::cerr << libphysica::Formatted_String("Error", "Red", true) << " in obscura::Atom::Atom(std::string): Element " << element_name << " not recognized." << std::endl;
		std::exit(EXIT_FAILURE);
	}
}

double Atom::Lowest_Binding_Energy() const
{
	double binding_energy_min = 1.0;
	for(unsigned int i = 0; i < electrons.size(); i++)
	{
		if(std::fabs(electrons[i].binding_energy) < binding_energy_min)
			binding_energy_min = std::fabs(electrons[i].binding_energy);
	}
	return binding_energy_min;
}

Atomic_Electron Atom::Electron(unsigned int n, unsigned int l)
{
	for(unsigned int i = 0; i < electrons.size(); i++)
	{
		if(electrons[i].n == n && electrons[i].l == l)
			return electrons[i];
	}
	std::cerr << libphysica::Formatted_String("Error", "Red", true) << " in obscura::Atom::Electron(): (n,l) = (" << n << "," << l << ") of " << nucleus.name << " does not exist." << std::endl;
	std::exit(EXIT_FAILURE);
}

void Atom::Print_Summary(unsigned int MPI_rank) const
{
	if(MPI_rank == 0)
	{
		std::cout << "Atom:\t\t" << nucleus.name << std::endl
				  << "Z:\t\t" << nucleus.Z << std::endl
				  << "Mass [GeV]:\t" << nucleus.Average_Nuclear_Mass() << std::endl
				  << "W [eV]:\t\t" << In_Units(W, eV) << std::endl
				  << "Electron orbitals:" << std::endl
				  << "\tOrbital\tE_Binding [eV]\tSecondary electrons" << std::endl;
		for(unsigned int i = 0; i < electrons.size(); i++)
		{
			std::cout << "\t";
			electrons[i].Print_Summary();
		}
	}
}

}	// namespace obscura
