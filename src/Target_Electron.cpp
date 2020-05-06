#include "Target_Electron.hpp"

#include <cmath>
#include <fstream>

//Headers from libphysica library
#include "Natural_Units.hpp"
#include "Utilities.hpp"


//1. Kinematic functions
	double vMinimal_Electrons(double q,double Delta_E, double mDM)
	{
		return (Delta_E / q + q / 2.0 / mDM);
	}

//2. Semiconductor crystal target
	Semiconductor::Semiconductor(std::string target)
	: name(target) , dE(0.1*eV) , dq(0.02*aEM*mElectron) 
	{
		double prefactor;
		if(name == "Si")
		{
			energy_gap = 1.11*eV;
			epsilon = 3.6*eV;
			M_cell = 2.0*28.08*mNucleon;
			prefactor = 2.0*eV;
			Q_max = std::floor( (50*eV - energy_gap + epsilon)/epsilon);
		}
		else if(name == "Ge")
		{
			energy_gap = 0.67*eV;
			epsilon = 2.9*eV;
			M_cell = 2.0*72.6*mNucleon;
			prefactor = 1.8*eV;
			Q_max = std::floor( (50*eV - energy_gap + epsilon)/epsilon);
		}
		else
		{
			std::cerr <<"Error in Semiconductor::Semiconductor(): target " <<target <<" not recognized."<<std::endl;
			std::exit(EXIT_FAILURE);
		}
		//Import the form factor
			std::ifstream f;
			f.open("../data/Semiconductors/C."+target+"137.dat");
			//Prefactors:
			double wk = 2.0/137.0;
			for(int Ei = 1; Ei <= 500; Ei++)
			{
				for(int qi=1;qi<=900;qi++)
				{
					f >> Crystal_Form_Factor[qi-1][Ei-1];
					Crystal_Form_Factor[qi-1][Ei-1] *= prefactor*qi/dE*wk/4.0;
				}
			}
			f.close();
	}

//3. Bound electrons in isolated atoms
	std::string s_names[5] = {"s","p","d","f","g"};

	Atomic_Electron::Atomic_Electron(std::string element,double A, int N, int L,double Ebinding, double kMin,double kMax, double qMin, double qMax, double w, unsigned int neSecondary)
	: k_min(kMin), k_max(kMax), q_min(qMin), q_max(qMax), n(N), l(L), binding_energy(Ebinding), nucleus_mass(A * mNucleon), W(w), number_of_secondary_electrons(neSecondary)
	{
		name = element + "_" + std::to_string(n) + s_names[l];
		//Import the table.
		Form_Factor_Tables = Import_Table("../data/Form_Factors_Ionization/"+name+".txt");
		Nk = Form_Factor_Tables.size();
		Nq = Form_Factor_Tables[0].size();
		k_Grid = Log_Space(k_min, k_max, Nk);
		q_Grid = Log_Space(q_min, q_max, Nq);
		dlogk = log10(k_max/k_min) / (Nk-1.0);
		dlogq = log10(q_max/q_min) / (Nq-1.0);
	}

	double Atomic_Electron::Ionization_Form_Factor(double q, double E) const
	{
		double k = sqrt(2.0*mElectron*E);
		int ki = std::floor(log10(k/k_min) / dlogk);
		int qi = std::floor(log10(q/q_min) / dlogq);
		if(ki == (int) Nk-1)		ki--;
		else if(ki == -1)	ki++;
		
		//Bilinear interpolation between four points: (Source: https://en.wikipedia.org/wiki/Bilinear_interpolation)
		double x = k;
		double y = q;
		double x1 = k_Grid[ki];
		double x2 = k_Grid[ki+1];
		double y1 = q_Grid[qi];
		double y2 = q_Grid[qi+1];
		double f11 = Form_Factor_Tables[ki][qi];
		double f12 = Form_Factor_Tables[ki][qi+1];
		double f21 = Form_Factor_Tables[ki+1][qi];
		double f22 = Form_Factor_Tables[ki+1][qi+1];
		
		double f = 1.0 / (x2-x1) / (y2 - y1) * ( f11*(x2-x)*(y2-y) + f21*(x-x1)*(y2-y) + f12*(x2-x)*(y-y1) + f22*(x-x1)*(y-y1));
		return f;
	}

	void Atomic_Electron::Print_Summary(unsigned int MPI_rank) const
	{
		if(MPI_rank == 0)
		{
			std::cout <<name<<"\t"<<In_Units(binding_energy,eV)<<"\t\t"<<number_of_secondary_electrons<<std::endl;
		}
	}		


	Atom::Atom(std::string element_name, int z, double a, double w, std::vector<Atomic_Electron> shells)
	: name(element_name), Z(z), A(a), mass(A * mNucleon), W(w), electrons(shells)
	{
	}

	double Atom::Lowest_Binding_Energy() const
	{
		double binding_energy_min = 1.0;
		for(unsigned int i = 0; i< electrons.size(); i++)
		{
			if( std::fabs(electrons[i].binding_energy) < binding_energy_min) binding_energy_min = std::fabs(electrons[i].binding_energy);
		}
		return binding_energy_min;
	}

	Atomic_Electron Atom::Electron(unsigned int n, unsigned int l)
	{
		for(unsigned int i = 0; i<electrons.size(); i++)
		{
			if(electrons[i].n == n && electrons[i].l == l) return electrons[i];
		}
		std::cerr <<"Error in Atom::Electron(): (n,l) = ("<<n<<","<<l<<") of " <<name<<" does not exist."<<std::endl;
		std::exit(EXIT_FAILURE);
	}

	void Atom::Print_Summary(unsigned int MPI_rank) const
	{
		if(MPI_rank == 0)
		{
			std::cout	<<"Atom:\t\t" <<name<<std::endl
						<<"(Z,A):\t\t(" <<Z<<","<<A<<")"<<std::endl
						<<"Mass [GeV]:\t"	<<mass <<std::endl
						<<"W [eV]:\t\t"	<<In_Units(W,eV)<<std::endl
						<<"Electron orbitals:"<<std::endl
						<<"\tOrbital\tE_Binding [eV]\tSecondary electrons"<<std::endl;
			for(unsigned int i = 0; i < electrons.size(); i++)
			{
				std::cout<<"\t";
				electrons[i].Print_Summary();
			}
		}
	}

	Atom Import_Ionization_Form_Factors(std::string element)
	{
		if(element == "Xe" || element == "Xenon")
		{
			//Xenon:
			int Z = 54;
			double A = 131.0;
			double W = 13.8*eV;
			double q_min = 1.0 * keV;
			double q_max = 1000.0 * keV;
			double k_min = 0.1 * keV;
			double k_max = 100.0 * keV;
			// Atomic_Electron Xe_1s("Xe", A, 1, 0, 33317.6*eV, k_min, k_max, q_min, q_max);
			// Atomic_Electron Xe_2s("Xe", A, 2, 0, 5149.21*eV, k_min, k_max, q_min, q_max);
			// Atomic_Electron Xe_2p("Xe", A, 2, 1, 4837.71*eV, k_min, k_max, q_min, q_max);
			// Atomic_Electron Xe_3s("Xe", A, 3, 0, 1093.24*eV, k_min, k_max, q_min, q_max);
			// Atomic_Electron Xe_3p("Xe", A, 3, 1, 958.43*eV, k_min, k_max, q_min, q_max);
			// Atomic_Electron Xe_3d("Xe", A, 3, 2, 710.73*eV, k_min, k_max, q_min, q_max);
			Atomic_Electron Xe_4s("Xe", A, 4, 0, 213.781*eV, k_min, k_max, q_min, q_max, W, 3);
			Atomic_Electron Xe_4p("Xe", A, 4, 1, 163.495*eV, k_min, k_max, q_min, q_max, W, 6);
			Atomic_Electron Xe_4d("Xe", A, 4, 2, 75.5897*eV, k_min, k_max, q_min, q_max, W, 4);
			Atomic_Electron Xe_5s("Xe", A, 5, 0, 25.6986*eV, k_min, k_max, q_min, q_max, W, 0);
			Atomic_Electron Xe_5p("Xe", A, 5, 1, 12.4433*eV, k_min, k_max, q_min, q_max, W, 0);
			return Atom("Xenon",Z,A,W,{Xe_5p,Xe_5s,Xe_4d,Xe_4p,Xe_4s});
		}
		else if (element == "Ar" || element == "Argon")
		{
			//Argon:
			int Z = 18;
			double A = 40.0;
			double W = 19.6*eV;
			double q_min = 1.0 * keV;
			double q_max = 1000.0 * keV;
			double k_min = 0.1 * keV;
			double k_max = 100.0 * keV;
			Atomic_Electron Ar_1s("Ar", A, 1, 0, 3227.55*eV, k_min, k_max, q_min, q_max, W, 0);
			Atomic_Electron Ar_2s("Ar", A, 2, 0, 335.303*eV, k_min, k_max, q_min, q_max, W, 0);
			Atomic_Electron Ar_2p("Ar", A, 2, 1, 260.453*eV, k_min, k_max, q_min, q_max, W, 0);
			Atomic_Electron Ar_3s("Ar", A, 3, 0, 34.7585*eV, k_min, k_max, q_min, q_max, W, 0);
			Atomic_Electron Ar_3p("Ar", A, 3, 1, 16.0824*eV, k_min, k_max, q_min, q_max, W, 0);
			return Atom("Argon",Z,A,W,{Ar_3p,Ar_3s,Ar_2p,Ar_2s,Ar_1s});
		}
		else
		{
			std::cerr <<"Error in Import_Ionization_Form_Factors(std::string): Element " <<element <<" not recognized."<<std::endl;
			std::exit(EXIT_FAILURE);
		}
	}
