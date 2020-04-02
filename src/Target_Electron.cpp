#include "Target_Electron.hpp"

#include <fstream>

//Headers from libphys library
#include "Natural_Units.hpp"

//1. Semiconductor crystal target
	//Semiconductor target
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
		}
		else if(name == "Ge")
		{
			energy_gap = 0.67*eV;
			epsilon = 2.9*eV;
			M_cell = 2.0*72.6*mNucleon;
			prefactor = 1.8*eV;
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
			for(int Ei=1;Ei<=500;Ei++)
			{
				for(int qi=1;qi<=900;qi++)
				{
					
					f >> Crystal_Form_Factor[qi-1][Ei-1];
					Crystal_Form_Factor[qi-1][Ei-1] *= prefactor*qi/dE*wk/4.0;
				}
			}
			f.close();
	}