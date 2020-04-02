#ifndef __Target_Electron_hpp_
#define __Target_Electron_hpp_

#include <string>

//1. Semiconductor crystal target
	struct Semiconductor 
	{
		std::string name;
		double dE,dq;
		double M_cell;
		double energy_gap, epsilon;
		double Crystal_Form_Factor[900][500];
		Semiconductor(std::string target);
	};

#endif