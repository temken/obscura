#include "Target_Electron.hpp"



// //2. Semiconductors
// 	//Crystal target
// 	Crystal::Crystal(std::string target)
// 	: name(target) , dE(0.1*eV) , dq(0.02*aEM*mElectron) 
// 	{
// 		double prefactor;
// 		if(name == "Si")
// 		{
// 			E_Gap = 1.11*eV;
// 			eps = 3.6*eV;
// 			Mcell = 2.0*28.08*mNucleon;
// 			prefactor = 2.0*eV;
// 		}
// 		else if(name == "Ge")
// 		{
// 			E_Gap = 0.67*eV;
// 			eps = 2.9*eV;
// 			Mcell = 2.0*72.6*mNucleon;
// 			prefactor = 1.8*eV;
// 		}
// 		else
// 		{
// 			std::cerr <<"Error in Crystal::Crystal(): target " <<target <<" not recognized."<<std::endl;
// 			std::exit(EXIT_FAILURE);
// 		}
// 		//Import the form factor
// 			std::ifstream f;
// 			f.open("../data/Semiconductors/C."+target+"137.dat");
// 			//Prefactors:
// 			double wk = 2.0/137.0;
// 			for(int Ei=1;Ei<=500;Ei++)
// 			{
// 				for(int qi=1;qi<=900;qi++)
// 				{
					
// 					f >> FormFactor[qi-1][Ei-1];
// 					FormFactor[qi-1][Ei-1] *= prefactor*qi/dE*wk/4.0;
// 				}
// 			}
// 			f.close();
// 	}