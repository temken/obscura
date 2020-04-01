#ifndef __Configuration__hpp_
#define __Configuration__hpp_

#include <string>

#include "Direct_Detection.hpp"
#include "DM_Distribution.hpp"
#include "DM_Particle.hpp"

struct Configuration
{
	std::string ID;

	//Mass scan
	double mass_min, mass_max;
	unsigned int number_of_masses;

	DM_Particle *DM =	{ nullptr };
	DM_Distribution *DM_distr =	{ nullptr };
	DM_Detector *DM_detector =	{ nullptr };

	//Constructors
	Configuration();
	Configuration(std::string cfg_filename, int MPI_rank = 0);
};

#endif