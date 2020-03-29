#ifndef __Target_Nucleus_hpp_
#define __Target_Nucleus_hpp_

#include <vector>
#include <string>

struct Isotope
{
		unsigned int Z,A;
		double abundance;
		double spin;
		double sp,sn;
		
		std::string name;
		double mass;

		Isotope();
		Isotope(unsigned int z, unsigned int a,double abund=1.0,double Spin=0.0,double Sp=0.0,double Sn=0.0);

		double Thomas_Fermi_Radius() const;
		double Helm_Form_Factor(double q) const;
};

class Element
{
	private:
		std::vector<Isotope> isotopes;

	public:
		std::string name;

		Element();
		Element(std::vector<Isotope>& iso);
		Element(const Isotope& iso);

		unsigned int Number_of_Isotopes() const;

		void Add_Isotope(Isotope& isotope);

		Isotope& operator[](int i) 
		{	
			return isotopes[i];
		}
		const Isotope& operator[](int i) const 
		{
			return isotopes[i];
		}

		double Average_Nuclear_Mass() const;

		void Print_Summary() const;
};

extern void Import_Nuclear_Data();

extern std::vector<Element> Elements;
extern Element Get_Element(int Z);
extern Element Get_Element(std::string name);

#endif