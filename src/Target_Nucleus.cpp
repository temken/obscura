#include "obscura/Target_Nucleus.hpp"

#include <cmath>
#include <fstream>
#include <iostream>

#include "libphysica/Natural_Units.hpp"
#include "libphysica/Special_Functions.hpp"
#include "libphysica/Utilities.hpp"

namespace obscura
{

using namespace libphysica::natural_units;

// 1. Kinematic functions
double vMinimal_Nucleus(double ER, double mDM, double mNucleus)
{
	return sqrt(mNucleus * ER / 2.0 / pow(libphysica::Reduced_Mass(mDM, mNucleus), 2.0));
}
double Maximum_Nuclear_Recoil_Energy(double vDM, double mDM, double mNucleus)
{
	return 2.0 * vDM * vDM * pow(libphysica::Reduced_Mass(mDM, mNucleus), 2.0) / mNucleus;
}

// 2. Class for nuclear isotopes.
// Auxiliary list with element names
std::vector<std::string> Nucleus_Names = {"H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"};

Isotope::Isotope()
: Z(1), A(1), abundance(1.0), spin(0.5), sp(0.5), sn(0)
{
	name = "H-1";
	mass = mProton;
}

Isotope::Isotope(unsigned int z, unsigned int a, double abund, double Spin, double Sp, double Sn)
: Z(z), A(a), abundance(abund), spin(Spin), sp(Sp), sn(Sn)
{
	name = Nucleus_Names[Z - 1] + "-" + std::to_string(A);
	mass = (A == 1) ? mProton : A * mNucleon;
}

double Isotope::Thomas_Fermi_Radius() const
{
	return pow(9 * M_PI * M_PI / 2.0 / Z, 1.0 / 3.0) / 4.0 * Bohr_Radius;
}

double Isotope::Helm_Form_Factor(double q) const
{
	if(q < 1.0e-6 * MeV)
		return 1.0;
	double a  = 0.52 * fm;
	double c  = (1.23 * pow(A, 1.0 / 3.0) - 0.6) * fm;
	double s  = 0.9 * fm;
	double rn = sqrt(c * c + 7.0 / 3.0 * pow(M_PI * a, 2.0) - 5.0 * s * s);
	double qr = q * rn;
	return 3.0 * (sin(qr) / pow(qr, 3.0) - cos(qr) / pow(qr, 2.0)) * exp(-q * q * s * s / 2.0);
}

void Isotope::Print_Summary(unsigned int MPI_rank) const
{
	if(MPI_rank == 0)
		std::cout << name << "\t" << Z << "\t" << A << "\t" << libphysica::Round(100.0 * abundance) << "\t\t" << spin << "\t" << sp << "\t" << sn << std::endl;
}

// 3. Class for elements containing all isotopes occuring in nature
Nucleus::Nucleus()
{
	isotopes = {};
}

Nucleus::Nucleus(const std::vector<Isotope>& iso)
: Z(iso[0].Z), isotopes(iso)
{
	name = Nucleus_Names[Z - 1];
}

Nucleus::Nucleus(const Isotope& iso)
: isotopes({iso})
{
	name = Nucleus_Names[iso.Z - 1];
}

unsigned int Nucleus::Number_of_Isotopes() const
{
	return isotopes.size();
}

// void Nucleus::Add_Isotope(Isotope& isotope)
// {
// 	isotopes.push_back(isotope);
// }

double Nucleus::Average_Nuclear_Mass() const
{
	double average_mass = 0.0;
	for(unsigned int i = 0; i < Number_of_Isotopes(); i++)
		average_mass += isotopes[i].mass * isotopes[i].abundance;
	return average_mass;
}

void Nucleus::Print_Summary(unsigned int MPI_rank) const
{
	if(MPI_rank == 0)
	{
		double total = 0.0;
		std::cout << std::endl
				  << name << std::endl
				  << "Isotope\tZ\tA\tAbund.[%]\tSpin\t<sp>\t<sn>" << std::endl;
		std::cout << "------------------------------------------------------------" << std::endl;
		for(const auto& isotope : isotopes)
		{
			total += 100.0 * isotope.abundance;
			isotope.Print_Summary(MPI_rank);
		}
		std::cout << "Total:\t\t" << libphysica::Round(Average_Nuclear_Mass() / mNucleon) << "\t" << total << std::endl
				  << std::endl;
	}
}

Isotope Nucleus::Get_Isotope(unsigned int A) const
{
	for(unsigned int i = 0; i < Number_of_Isotopes(); i++)
		if(isotopes[i].A == A)
			return isotopes[i];
	std::cerr << libphysica::Formatted_String("Error", "Red", true) << " in obscura::Nucleus::Get_Isotope(): Isotope A=" << A << " not existent for " << name << "." << std::endl;
	std::exit(EXIT_FAILURE);
}

// 4. Nuclear data
// Import the nuclear data from a file
std::vector<Nucleus> all_nuclei;
std::vector<Nucleus> Import_Nuclear_Data()
{
	std::vector<Nucleus> nuclei = {};
	std::vector<Isotope> isotopes;
	std::string path = PROJECT_DIR "data/Nuclear_Data.txt";
	std::ifstream f;
	f.open(path);
	if(!f)
	{
		std::cerr << libphysica::Formatted_String("Error", "Red", true) << " in obscura::Import_Nuclear_Data(): Data file " << path << " not found." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	std::string name;
	int Z, A;
	double abund, spin, sp, sn;
	int Zold = 1;
	while(f >> name >> Z >> A >> abund >> spin >> sp >> sn)
	{
		if(Z > Zold)
		{
			nuclei.push_back(Nucleus(isotopes));
			isotopes.clear();
			Zold = Z;
		}
		isotopes.push_back(Isotope(Z, A, abund, spin, sp, sn));
	}
	nuclei.push_back(Nucleus(isotopes));
	f.close();
	return nuclei;
}

Isotope Get_Isotope(unsigned int Z, unsigned int A)
{
	return Get_Nucleus(Z).Get_Isotope(A);
}

Nucleus Get_Nucleus(unsigned int Z)
{
	if(Z < 1 || Z > 92)
	{
		std::cerr << libphysica::Formatted_String("Error", "Red", true) << " in obscura::Get_Nucleus(): Input Z=" << Z << " is not a value between 1 and 92." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	else if(all_nuclei.size() == 0)
		all_nuclei = Import_Nuclear_Data();

	return all_nuclei[Z - 1];
}

Nucleus Get_Nucleus(std::string name)
{
	for(int Z = 1; Z <= 92; Z++)
	{
		if(Get_Nucleus(Z).name == name)
			return Get_Nucleus(Z);
	}
	std::cerr << libphysica::Formatted_String("Error", "Red", true) << " in obscura::Get_Nucleus(): Nucleus " << name << " not recognized." << std::endl;
	std::exit(EXIT_FAILURE);
}

}	// namespace obscura
