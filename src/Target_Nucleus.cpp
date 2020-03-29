#include "Target_Nucleus.hpp"

#include <iostream>
#include <fstream>
#include <cmath>

//Headers from libphys library
#include "Numerics.hpp"
#include "Natural_Units.hpp"

//Auxiliary list with element names
	std::vector<std::string> ElementNames={"H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg","Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr","Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br","Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd","Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La","Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er","Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au","Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th","Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md","No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn","Nh", "Fl", "Mc", "Lv", "Ts", "Og"};

//Class for isotopes
	Isotope::Isotope()
	: Z(1), A(1), abundance(1.0), spin(0.5), sp(0.5), sn(0)
	{
		name="H-1";
		mass = mProton;
	}

	Isotope::Isotope(unsigned int z, unsigned int a,double abund,double Spin,double Sp,double Sn)
	: Z(z), A(a), abundance(abund), spin(Spin), sp(Sp), sn(Sn)
	{
		name=ElementNames[Z-1]+"-"+std::to_string(A);
		mass= (A==1)? mProton : A*mNucleon;
	}

	double Isotope::Thomas_Fermi_Radius() const
	{
		return pow(9*M_PI*M_PI/2.0/Z,1.0/3.0)/4.0*Bohr_Radius;
	}
	
	double Isotope::Helm_Form_Factor(double q) const
	{
		double a = 0.52*fm;
		double c = (1.23 * pow(A,1.0/3.0) - 0.6)*fm;
		double s = 0.9*fm;
		double rn = sqrt(c*c + 7.0/3.0*pow(M_PI*a,2.0) - 5.0*s*s);
		double qr = q*rn;
		return 3.0 * ( sin(qr) / pow(qr,3.0) - cos(qr) / pow(qr,2.0) ) * exp(-q*q*s*s/2.0);
	}

	

//Class for elements containing all isotopes occuring in nature
	Element::Element()
	{
		isotopes = {};
	}

	Element::Element(std::vector<Isotope>& iso)
	: isotopes(iso)
	{
		name=ElementNames[iso[0].Z-1];
	}

	Element::Element(const Isotope& iso)
	: isotopes({iso})
	{
		name=ElementNames[iso.Z-1];
	}

	unsigned int Element::Number_of_Isotopes() const
	{
		return isotopes.size();
	}

	void Element::Add_Isotope(Isotope& isotope)
	{
		isotopes.push_back(isotope);
	}

	double Element::Average_Nuclear_Mass() const
	{
		double average_mass=0.0;
		for(unsigned int i=0;i<Number_of_Isotopes();i++) average_mass += isotopes[i].mass*isotopes[i].abundance;
		return average_mass;
	}


	void Element::Print_Summary() const
	{
		double total=0.0;
		std::cout <<std::endl<<name<<std::endl<<"Isotope\tZ\tA\tAbund.[%]\tSpin\t<sp>\t<sn>"<<std::endl;
		std::cout <<"------------------------------------------------------------"<<std::endl;
		for(unsigned int i=0;i<Number_of_Isotopes();i++)
		{
			total+=100.0*isotopes[i].abundance;
			std::cout <<isotopes[i].name <<"\t"<<isotopes[i].Z<<"\t"<<isotopes[i].A<<"\t"<<100.0*isotopes[i].abundance<<"\t\t"<<isotopes[i].spin<<"\t"<<isotopes[i].sp<<"\t"<<isotopes[i].sn<<std::endl;
		}
		std::cout <<"Total:\t\t"<<Average_Nuclear_Mass()/mNucleon<<"\t"<<total<<std::endl<<std::endl;
	}

	//Import the nuclear data from a file
	std::vector<Element> Elements;
	void Import_Nuclear_Data()
	{
		Elements.clear();
		std::vector<Isotope> isotopes;
		std::ifstream f;
		f.open("../data/Nuclear_Data.txt");

		std::string name;
		int Z,A;
		double abund,spin,sp,sn;
		int Zold=1;
		while(f >>name >>Z >>A >>abund >>spin >>sp >>sn)
		{
			if(Z>Zold)
			{
				Elements.push_back(Element(isotopes));
				isotopes.clear();
				Zold=Z;
			}
			isotopes.push_back( Isotope(Z, A,abund,spin,sp, sn) );
		}
		Elements.push_back(Element(isotopes));
		f.close();
	}

	Element Get_Element(int Z)
	{
		if(Z<1||Z>92)
		{
			std::cerr <<"Error in Get_Element(): Input Z="<<Z <<" is not a value between 1 and 92."<<std::endl;
			std::exit(EXIT_FAILURE);
		}
		else if(Elements.size() == 0)
		{
			std::cerr <<"Warning in Get_Element(): Nuclear data not imported yet. Import_Nuclear_Data() is performed now."<<std::endl;
			Import_Nuclear_Data();
		}
		return Elements[Z-1];
	}

	Element Get_Element(std::string name)
	{
		for(int Z = 1; Z<=92 ; Z++)
		{
			if(Get_Element(Z).name == name) return Get_Element(Z);
		}
		std::cerr<<"Error in Get_Element(): Element " <<name <<" not recognized."<<std::endl;
		std::exit(EXIT_FAILURE);
	}
