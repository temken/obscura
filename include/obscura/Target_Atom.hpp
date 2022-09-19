#ifndef __Target_Atom_hpp_
#define __Target_Atom_hpp_

#include <string>
#include <vector>

#include "libphysica/Numerics.hpp"

#include "Target_Nucleus.hpp"
#include "version.hpp"

namespace obscura
{

// 1. Kinematic functions
extern double vMinimal_Electrons(double q, double Delta_E, double mDM);

// 2. Bound electrons in isolated atoms
struct Atomic_Electron
{
	// Ionization form factor tables
	double k_min, k_max, q_min, q_max;
	double dlogk, dlogq;
	unsigned int Nk, Nq;
	std::vector<double> k_Grid = {};
	std::vector<double> q_Grid = {};
	std::vector<libphysica::Interpolation_2D> atomic_response_interpolations;

	unsigned int n, l;
	std::string name;
	double binding_energy;
	unsigned int number_of_secondary_electrons;

	Atomic_Electron(std::string element, int N, int L, double Ebinding, double kMin, double kMax, double qMin, double qMax, unsigned int neSecondary = 0);

	double Atomic_Response_Function(int response, double q, double E);
	// Squared ionization form factor.
	double Ionization_Form_Factor(double q, double E);

	void Print_Summary(unsigned int MPI_rank = 0) const;
};

struct Atom
{
	Nucleus nucleus;
	std::vector<Atomic_Electron> electrons;
	double W;

	// Constructor
	Atom(int z, double w, std::vector<Atomic_Electron> shells = {});
	Atom(const std::string& element_name);

	double Lowest_Binding_Energy() const;

	Atomic_Electron Electron(unsigned int n, unsigned int l);

	// Overloading brackets
	Atomic_Electron& operator[](int i)
	{
		return electrons[i];
	}

	const Atomic_Electron& operator[](int i) const
	{
		return electrons[i];
	}

	void Print_Summary(unsigned int MPI_rank = 0) const;
};

}	// namespace obscura

#endif