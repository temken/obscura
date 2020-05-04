#ifndef __Target_Electron_hpp_
#define __Target_Electron_hpp_

#include <string>
#include <vector>

//1. Kinematic functions
	extern double vMinimal_Electrons(double q,double Delta_E, double mDM);
	
//2. Semiconductor crystal target
	struct Semiconductor 
	{
		std::string name;
		double dE,dq;
		double M_cell;
		double energy_gap, epsilon;
		unsigned int Q_max;
		double Crystal_Form_Factor[900][500];
		Semiconductor(std::string target);
	};

//3. Bound electrons in isolated atoms
	struct Atomic_Electron
	{
			// Ionization form factor tables
			double k_min,k_max,q_min,q_max;
			double dlogk, dlogq;
			unsigned int Nk, Nq;
			std::vector<double> k_Grid={};
			std::vector<double> q_Grid={};
			std::vector< std::vector<double>> Form_Factor_Tables;

			unsigned int n,l;
			std::string name;
			double binding_energy;
			double nucleus_mass;
			double W;
			unsigned int number_of_secondary_electrons;

			Atomic_Electron(std::string element,double A, int N, int L, double Ebinding, double kMin,double kMax, double qMin, double qMax, double w, unsigned int neSecondary = 0);

			//Squared ionization form factor.
			double Ionization_Form_Factor(double q, double E) const;

			void Print_Summary(unsigned int MPI_rank = 0) const;		
	};

	struct Atom
	{
		std::string name;

		int Z;
		double A;
		
		double mass;

		double W;

		std::vector<Atomic_Electron> electrons;

		//Constructor
		Atom(std::string element_name, int z, double a, double w, std::vector<Atomic_Electron> shells = {});

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

	extern Atom Import_Ionization_Form_Factors(std::string element);


#endif